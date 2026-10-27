"""Figure F: negative weights (C12 by channel; pieces; nuclei)."""
from common import *
from matplotlib.lines import Line2D

MINCNT = 200
Q_EDGES = np.arange(0, 1.2001e3, 40.0)
RATIO_EDGES = np.arange(0.9, 2.0001, 0.02)
HR_EDGES = np.arange(-4, 8.0001, 0.1)


def binned_neg(x, ws, edges, mask=None):
    """per slot: n events, n negative, sum|w| of negatives, sum w.   ws = list of 3 weight arrays"""
    b = Binner(x, edges, mask)
    n = b.count()
    out = {"n": n, "nneg": [], "sneg": [], "ssum": []}
    for w in ws:
        neg = w < 0
        out["nneg"].append(b.count(neg))
        out["sneg"].append(b.sums(np.where(neg, -w, 0.0))[0])
        out["ssum"].append(b.sums(w)[0])
    return out


R = {}   # (nuc, beam) -> reduction
for nuc in NUCS:
    for beam in BEAMS:
        s = load_sample(nuc, beam)
        if s is None:
            continue
        ws = [s.w(k) for k in range(3)]
        wt = s.w()
        red = {"s": {"A": s.A, "Z": s.Z, "N": s.N}}
        eb = s.ebins()
        red["E"] = (eb, binned_neg(s.E, ws, eb))
        if nuc == "C12":
            red["q0"] = (Q_EDGES, binned_neg(s.q0, ws, Q_EDGES))
            red["qz"] = (Q_EDGES, binned_neg(s.qz, ws, Q_EDGES))
        if nuc in ("Ar40", "Pb208"):
            rmax = np.ceil(s.a("r").max() / 0.3) * 0.3
            red_r = np.arange(0, rmax + 1e-9, 0.3)
            red["r"] = (red_r, binned_neg(s.a("r"), ws, red_r))
            with np.errstate(divide="ignore", invalid="ignore"):
                kr = s.a("kF_n") / s.a("kF_p")
            ok = np.isfinite(kr)
            red["kr"] = (RATIO_EDGES, binned_neg(np.where(ok, kr, -1), ws, RATIO_EDGES))
            red["kr_stats"] = {"min": float(np.nanmin(kr[ok])), "p99.5": float(np.nanpercentile(kr[ok], 99.5)), "max": float(np.nanmax(kr[ok]))}
        # scalars
        sig = sum(w.sum() for w in ws)
        red["scal"] = {
            "neg_event_frac_pooled": float(np.mean([(w < 0).mean() for w in ws])),
            "neg_event_frac_per_slot": [float((w < 0).mean()) for w in ws],
            "neg_event_frac_total_weight": float((wt < 0).mean()),
            "neg_xsec_frac_pooled": float(sum(np.where(w < 0, -w, 0).sum() for w in ws) / sig),
            "neg_xsec_frac_per_slot": [float(np.where(w < 0, -w, 0).sum() / w.sum()) for w in ws],
            "sigma_tot": float(sig * s.fac),
        }
        pieces = {}
        for p in ("rho", "pi", "int"):
            pieces[p] = float(sum(np.minimum(s.a(f"w{k}_{p}"), 0).sum() for k in range(3)) / sig)
        red["piece_neg_over_sigma"] = pieces
        red["piece_pos_over_sigma"] = {p: float(sum(np.maximum(s.a(f"w{k}_{p}"), 0).sum() for k in range(3)) / sig) for p in ("rho", "pi", "int")}
        # distribution of w_int/w for events (slot-level) with w<0
        h = np.zeros(len(HR_EDGES) - 1)
        under = over = tot_neg = 0
        for k in range(3):
            m = ws[k] < 0
            ratio = s.a(f"w{k}_int")[m] / ws[k][m]
            h += np.histogram(np.clip(ratio, HR_EDGES[0], HR_EDGES[-1] - 1e-9), HR_EDGES)[0]
            under += (ratio < HR_EDGES[0]).sum(); over += (ratio > HR_EDGES[-1]).sum(); tot_neg += m.sum()
            red.setdefault("int_ratio_stats", []).append({"median": float(np.median(ratio)), "frac_int_negative": float((s.a(f"w{k}_int")[m] < 0).mean()),
                                                           "frac_ratio_gt_1": float((ratio > 1).mean())})
        red["hist_int_ratio"] = (h, int(under), int(over), int(tot_neg))
        R[(nuc, beam)] = red
        del ws, wt, s

extra = {f"{n}_{b}": r["scal"] | {"piece_neg_over_sigma": r["piece_neg_over_sigma"], "piece_pos_over_sigma": r["piece_pos_over_sigma"]} for (n, b), r in R.items()}
extra_i = {f"{n}_{b}": r["int_ratio_stats"] for (n, b), r in R.items()}
update_summary({"F_neg_all": extra, "F_neg_int_ratio_stats": extra_i})
miss = missing_list()


def line(ax, xc, num, den, mask, **kw):
    with np.errstate(divide="ignore", invalid="ignore"):
        y = np.where(mask, num / den, np.nan)
    ax.plot(xc, y, **kw)


# ----------------------------------------------------------------------------- F_neg_C12
fig, axs = plt.subplots(2, 2, figsize=(14, 10))
panels = [("E", "n_frac", r"$E_\nu$ [GeV]", "fraction of events with $w<0$", 1e-3),
          ("E", "x_frac", r"$E_\nu$ [GeV]", r"$\sum_{w<0}|w| \,/\, \sum w$", 1e-3),
          ("q0", "x_frac", r"$q_0$ [GeV]", r"$\sum_{w<0}|w| \,/\, \sum w$", 1e-3),
          ("qz", "x_frac", r"$|\vec q|$ [GeV]", r"$\sum_{w<0}|w| \,/\, \sum w$", 1e-3)]
for ax, (var, kind, xl, yl, sc), lab in zip(axs.flat, panels, "abcd"):
    drawn = False
    for beam in BEAMS:
        if ("C12", beam) not in R or var not in R[("C12", beam)]:
            continue
        edges, d = R[("C12", beam)][var]
        xc = centres(edges) * sc
        ok = d["n"] >= MINCNT
        for k in range(3):
            if kind == "n_frac":
                line(ax, xc, d["nneg"][k], d["n"], ok, color=SLOT_COL[k], ls=BEAM_LS[beam], marker=BEAM_MARK[beam], ms=3)
            else:
                line(ax, xc, d["sneg"][k], d["ssum"][k], ok & (d["ssum"][k] > 0), color=SLOT_COL[k], ls=BEAM_LS[beam], marker=BEAM_MARK[beam], ms=3)
        drawn = True
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    ax.set_title(f"({lab}) " + {"n_frac": "negative-weight event fraction", "x_frac": "negative fraction of the cross section"}[kind] + {"E": " vs $E_\\nu$", "q0": " vs $q_0$", "qz": " vs $|\\vec q|$"}[var], fontsize=10.5)
    if not drawn:
        no_data(ax, "C12 samples missing")
    elif kind == "x_frac":
        ax.set_yscale("log"); ax.set_ylim(1e-4, None)
h = [Line2D([], [], color=SLOT_COL[k], lw=2, label=mirror_label(k)) for k in range(3)]
h += [Line2D([], [], color="k", ls="-", marker="o", ms=3, label=r"$\nu$"), Line2D([], [], color="k", ls="--", marker="s", ms=3, label=r"$\bar\nu$")]
fig.legend(handles=h, fontsize=9.5, loc="lower center", ncol=5, bbox_to_anchor=(0.5, 0.02))
fig.suptitle(r"F. Negative weights in $^{12}$C, per final-state channel (bins with fewer than %d events omitted)" % MINCNT, fontsize=12)
stamp_missing(fig, [m for m in miss if m.startswith("C12")])
fig.tight_layout(rect=(0, 0.06, 1, 0.96))
save_fig(fig, "F_neg_C12", "12C negative weights per channel: fraction of events with w<0 vs E_nu (a) and negative fraction of the cross section sum_{w<0}|w|/sum(w) vs E_nu (b), q0 (c) and |q| (d); nu solid, nubar dashed, mirror channels (pp (nu) / nn (nubar), np (nu) / pn (nubar), pn (nu) / np (nubar)) in the same colour.")

# ----------------------------------------------------------------------------- F_neg_pieces
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(17, 6.2), gridspec_kw={"width_ratios": [1.6, 1, 1]})
keys = [(n, b) for n in NUCS for b in BEAMS if (n, b) in R]
if keys:
    x = np.arange(len(keys))
    bot = np.zeros(len(keys))
    for p in ("rho", "pi", "int"):
        v = np.array([R[k]["piece_neg_over_sigma"][p] for k in keys])
        ax1.bar(x, v, 0.75, bottom=bot, color=PIECE_COL[p], label=PIECE_LABEL[p], edgecolor="white", lw=0.4)
        bot += v
    ax1.set_xticks(x); ax1.set_xticklabels([f"{n}\n{'nu' if b == 'nu' else 'nubar'}" for n, b in keys], fontsize=9)
    ax1.axhline(0, color="k", lw=0.8)
    ax1.set_ylabel(r"$\sum \min(w_{piece},0)\,/\,\sum w$  (summed over all channels)")
    ax1.set_title("(a) negative contribution of each piece, relative to total $\\sigma$ (stacked)")
    ax1.set_ylim(bottom=ax1.get_ylim()[0] * 1.3)
    ax1.legend(loc="lower right", fontsize=9)
else:
    no_data(ax1)
for ax, beam, lab in [(ax2, "nu", "b"), (ax3, "nubar", "c")]:
    drawn = False
    xc = centres(HR_EDGES)
    for nuc in NUCS:
        if (nuc, beam) in R:
            h, under, over, tot = R[(nuc, beam)]["hist_int_ratio"]
            ax.step(xc, h / max(tot, 1) / np.diff(HR_EDGES), where="mid", color=NUC_COL[nuc], label=f"{nuc}  (<{HR_EDGES[0]:.0f}: {100*under/max(tot,1):.1f}%, >{HR_EDGES[-1]:.0f}: {100*over/max(tot,1):.1f}%)")
            drawn = True
    ax.axvline(1, color="k", lw=0.6, ls=":"); ax.axvline(0, color="k", lw=0.6, ls=":")
    ax.set_yscale("log")
    ax.set_xlabel(r"$w_{int}/w$   per channel, for events with $w<0$")
    ax.set_ylabel("probability density")
    ax.set_title(f"({lab}) interference share of negative weights, {'$\\nu$' if beam == 'nu' else '$\\\\bar\\\\nu$'}".replace("\\\\", "\\"), fontsize=10.5)
    if drawn:
        ax.legend(fontsize=9, title="out-of-range share", title_fontsize=9)
    else:
        no_data(ax)
stamp_missing(fig, miss)
fig.suptitle("F. Which piece makes the weights negative?  (w = w_rho + w_pi + w_int;  w_int/w > 1 means the interference term alone drives w below 0)", fontsize=11.5)
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
save_fig(fig, "F_neg_pieces", "Negative contributions by piece (rho-transverse, pion-longitudinal, interference) relative to the total cross section for every sample (a), and the distribution of w_int/w for events with w<0 (b nu, c nubar).")

# ----------------------------------------------------------------------------- F_neg_nuclei
fig, axs = plt.subplots(3, 2, figsize=(13.5, 14))
asym = lambda n: (N_OF[n] - Z_OF[n]) / A_OF[n]
for j, (key, yl) in enumerate([("neg_event_frac_pooled", "fraction of (event, channel) pairs with $w<0$"), ("neg_xsec_frac_pooled", r"$\sum_{w<0}|w|/\sum w$  (all channels)")]):
    ax = axs[0, j]
    for (n, b), r in R.items():
        dx = -0.003 if b == "nu" else 0.003
        ax.scatter(asym(n) + dx, r["scal"][key], s=70, color=NUC_COL[n], marker=BEAM_MARK[b], edgecolor="k", lw=0.5, zorder=3)
        if b == "nu":
            ax.annotate(n, (asym(n) + dx, r["scal"][key]), textcoords="offset points", xytext=(-6, 5), fontsize=9, ha="right")
    ax.set_xlabel(r"$(N-Z)/A$"); ax.set_ylabel(yl)
    ax.set_title("(%s) %s vs neutron excess" % ("ab"[j], "negative event fraction" if j == 0 else "negative cross-section fraction"), fontsize=10.5)
    if not R:
        no_data(ax)
hh = [Line2D([], [], color="gray", marker="o", ls="", mec="k", label=r"$\nu$"), Line2D([], [], color="gray", marker="s", ls="", mec="k", label=r"$\bar\nu$")]
axs[0, 0].legend(handles=hh, loc="upper center", ncol=2, fontsize=9)
for row, (var, xl) in enumerate([("r", "$r$ [fm]"), ("kr", r"local $k_{F,n}/k_{F,p}$")], start=1):
    for j, kind in enumerate(["n_frac", "x_frac"]):
        ax = axs[row, j]
        drawn = False
        for nuc in ("Ar40", "Pb208"):
            for beam in BEAMS:
                if (nuc, beam) not in R:
                    continue
                edges, d = R[(nuc, beam)][var]
                xc = centres(edges)
                ok = d["n"] >= MINCNT
                nneg = sum(d["nneg"]); n = 3 * d["n"]
                sneg = sum(d["sneg"]); ss = sum(d["ssum"])
                if kind == "n_frac":
                    line(ax, xc, nneg, n, ok, color=NUC_COL[nuc], ls=BEAM_LS[beam], marker=BEAM_MARK[beam], ms=3, label=f"{nuc} {beam}")
                else:
                    line(ax, xc, sneg, ss, ok & (ss > 0), color=NUC_COL[nuc], ls=BEAM_LS[beam], marker=BEAM_MARK[beam], ms=3, label=f"{nuc} {beam}")
                drawn = True
        ax.set_xlabel(xl)
        ax.set_ylabel("fraction of (event, channel) pairs with $w<0$" if kind == "n_frac" else r"$\sum_{w<0}|w|/\sum w$  (all channels)")
        ax.set_title(f"({'cdef'[(row - 1) * 2 + j]}) vs {'radius' if var == 'r' else 'local kF ratio'}", fontsize=10.5)
        if drawn:
            ax.legend(fontsize=9)
        else:
            no_data(ax, "Ar40/Pb208 samples missing")
stamp_missing(fig, miss)
fig.suptitle("F. Negative weights across nuclei (pooled per channel: each event counted once per channel)", fontsize=12)
fig.tight_layout(rect=(0, 0.02, 1, 0.97))
save_fig(fig, "F_neg_nuclei", "Negative-weight fractions per nucleus and beam vs (N-Z)/A (a,b), and vs radius (c,d) and local kF_n/kF_p (e,f) for Ar40 and Pb208.")
finish("F")
