"""Figure B: hadron-tensor level mirror-symmetry unit tests (harness CSVs)."""
import pandas as pd
from common import *


def load(fn):
    p = HARNESS / fn
    if not p.exists():
        note_missing(f"harness/{fn}")
        return None
    return pd.read_csv(p, dtype={"comp": str})


def value(d):
    """physical value: Re part, except comp 12 where Im is physical (W3)."""
    is12 = (d["comp"] == "12").to_numpy()
    a = np.where(is12, d["a_im"], d["a_re"])
    b = np.where(is12, d["b_im"], d["b_re"])
    return a, b


def normalise(d, keys=("point", "param", "piece", "comp")):
    """returns a/max|a|, b/max|a|, with max over the three slots in each (point,param,piece,comp) group."""
    d = d.copy()
    d["a"], d["b"] = value(d)
    mx = d.groupby(list(keys))["a"].transform(lambda s: np.abs(s).max())
    d = d[mx > 0].copy()
    mx = mx[mx > 0]
    d["an"] = d["a"] / mx
    d["bn"] = d["b"] / mx
    d["dn"] = np.abs(d["bn"] - d["an"])
    return d


new = load("mirror_tests.csv")
old = load("mirror_tests_old.csv")
fig, axs = plt.subplots(2, 2, figsize=(13, 11.5))
extra = {}


def diag_panel(ax, d, title, grey=None):
    if grey is not None:
        ax.scatter(grey["an"], grey["bn"], s=3, c="#b0b0b0", alpha=0.35, lw=0, label=r"control: $\bar\nu$ on the un-mirrored gas (must fail)", rasterized=True)
    for k in range(3):
        m = d["slot"] == k
        ax.scatter(d.loc[m, "an"], d.loc[m, "bn"], s=4, c=SLOT_COL[k], alpha=0.5, lw=0,
                   label=f"{SLOT_FINAL['nu'][k]} ($\\nu$) ↔ {SLOT_FINAL['nubar'][k]} ($\\bar\\nu$)", rasterized=True)
    ax.plot([-1, 1], [-1, 1], "k-", lw=0.8)
    ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.05, 1.05)
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.legend(loc="upper left", markerscale=3, fontsize=9, frameon=True, framealpha=0.85, edgecolor="none")


# (a) T1_sym new ------------------------------------------------------------------
ax = axs[0, 0]
if new is not None:
    d = normalise(new[new["test"] == "T1_sym"])
    diag_panel(ax, d, r"(a) Joanna's branch: $\bar\nu$ vs $\nu$, N = Z matter")
    ax.set_xlabel(r"$\nu$: $W^{00}, W^{11}, W^{33}, \mathrm{Re}\,W^{03}, \mathrm{Im}\,W^{12}$ (normalised)"); ax.set_ylabel(r"$\bar\nu$, mirror channel: same component (normalised)")
    extra["T1_sym_new"] = {"max_abs_dev": float(d["dn"].max()), "median_abs_dev": float(d["dn"].median()),
                           "n_points": int(len(d)),
                           "max_dev_per_slot": [float(d.loc[d.slot == k, "dn"].max()) for k in range(3)]}
    ax.text(0.97, 0.05, f"max relative difference = {d['dn'].max():.1e}", transform=ax.transAxes, ha="right", fontsize=9.5)
    ax.text(0.03, 0.62, "each dot: one hadron-tensor component\n"
            r"($W^{00}$ charge, $W^{33}$ longitudinal current, $\mathrm{Re}\,W^{03}$ their" "\n"
            r"interference, $W^{11}$ transverse current, $\mathrm{Im}\,W^{12}$ V–A)" "\n"
            "for one channel and one piece (ρ, π, interference),\n"
            "at one random point (ω, |q|, pion and hole momenta),\n"
            "divided by the largest |value| of that component", transform=ax.transAxes, fontsize=8.5, va="top",
            bbox=dict(facecolor="white", edgecolor="0.8", alpha=0.9))
else:
    no_data(ax)
# (b) T1_sym_old --------------------------------------------------------------------
ax = axs[0, 1]
if old is not None:
    d = normalise(old[old["test"] == "T1_sym_old"])
    diag_panel(ax, d, "(b) Old code: same test")
    ax.set_xlabel(r"$\nu$: $W^{00}, W^{11}, W^{33}, \mathrm{Re}\,W^{03}, \mathrm{Im}\,W^{12}$ (normalised)"); ax.set_ylabel(r"$\bar\nu$, mirror channel: same component (normalised)")
    extra["T1_sym_old"] = {"max_abs_dev": float(d["dn"].max()), "median_abs_dev": float(d["dn"].median()),
                           "max_dev_per_slot": [float(d.loc[d.slot == k, "dn"].max()) for k in range(3)],
                           "median_dev_per_slot": [float(d.loc[d.slot == k, "dn"].median()) for k in range(3)],
                           "frac_off_diag_gt_1e-6_per_slot": [float((d.loc[d.slot == k, "dn"] > 1e-6).mean()) for k in range(3)]}
    nout = int(((d["an"].abs() > 1.05) | (d["bn"].abs() > 1.05)).sum())
    ax.text(0.97, 0.03, f"max relative difference = {d['dn'].max():.2f}\n{nout} of {len(d)} points outside the plotted range", transform=ax.transAxes, ha="right", fontsize=9.5,
            bbox=dict(fc="white", ec="none", alpha=0.8))
else:
    no_data(ax)
# (c) T2_asym + negative control -------------------------------------------------------
ax = axs[1, 0]
if new is not None:
    d = normalise(new[new["test"] == "T2_asym"])
    g = normalise(new[new["test"] == "NC_asym_unmirrored"])
    diag_panel(ax, d, r"(c) Joanna's branch, asymmetric matter (Ar40, Pb208 densities):" "\n" r"$\nu$ on ($k_{Fp},k_{Fn}$) vs $\bar\nu$ on ($k_{Fn},k_{Fp}$)", grey=g)
    ax.set_xlabel(r"$\nu$: $W^{00}, W^{11}, W^{33}, \mathrm{Re}\,W^{03}, \mathrm{Im}\,W^{12}$ (normalised)"); ax.set_ylabel(r"$\bar\nu$, mirror channel: same component (normalised)")
    extra["T2_asym"] = {"max_abs_dev": float(d["dn"].max()), "median_abs_dev": float(d["dn"].median()),
                        "param_kFn_over_kFp_range": [float(new[new.test == "T2_asym"].param.min()), float(new[new.test == "T2_asym"].param.max())]}
    extra["NC_asym_unmirrored"] = {"max_abs_dev": float(g["dn"].max()), "median_abs_dev": float(g["dn"].median()),
                                   "frac_dev_gt_1e-3": float((g["dn"] > 1e-3).mean())}
    ax.text(0.97, 0.03, f"mirrored: max relative difference {d['dn'].max():.1e}\ncontrol: median relative difference {g['dn'].median():.2f}, {100*(g['dn'] > 1e-3).mean():.0f}% of points > 1e-3", transform=ax.transAxes,
            ha="right", fontsize=9.5, bbox=dict(fc="white", ec="none", alpha=0.8))
else:
    no_data(ax)
# (d) T3 ------------------------------------------------------------------------------
ax = axs[1, 1]
if new is not None:
    from matplotlib.lines import Line2D
    FLOOR = 1e-19
    ep_all = []
    for test, ls, mk in [("T3_nu", "-", "o"), ("T3_nubar", "--", "s")]:
        t = new[new["test"] == test]
        if t.empty:
            continue
        d = normalise(t)
        for piece in ["rho", "pi", "int"]:
            dd = d[d.piece == piece]
            eps = np.sort(dd["param"].unique())
            med = np.array([dd.loc[dd.param == e, "dn"].median() for e in eps])
            p90 = np.array([dd.loc[dd.param == e, "dn"].quantile(0.9) for e in eps])
            frac = np.array([(dd.loc[dd.param == e, "dn"] > 1e-6).mean() for e in eps])
            ax.plot(eps, np.maximum(med, FLOOR), ls=ls, marker=mk, color=PIECE_COL[piece])
            ax.plot(eps, np.maximum(p90, FLOOR), ls=":", marker=mk, mfc="none", color=PIECE_COL[piece])
            extra.setdefault("T3", {})[f"{test}_{piece}"] = {"eps": eps, "median": med, "p90": p90,
                                                              "frac_dev_gt_1e-6": frac}
            ep_all += list(eps)
    if ep_all:
        e = np.array(sorted(set(ep_all)))
        ax.plot(e, e, "k--", lw=1)
    h = [Line2D([], [], color=PIECE_COL[p], lw=2, label=PIECE_LABEL[p]) for p in ["rho", "pi", "int"]]
    h += [Line2D([], [], color="k", ls="-", marker="o", label=r"$\nu$ (median: filled, solid)"),
          Line2D([], [], color="k", ls="--", marker="s", label=r"$\bar\nu$"),
          Line2D([], [], color="k", ls=":", marker="o", mfc="none", label="90th percentile (open, dotted)"),
          Line2D([], [], color="k", ls="--", lw=1, label=r"$\propto\epsilon$")]
    ax.legend(handles=h, loc="center", bbox_to_anchor=(0.5, 0.33), fontsize=9)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"relative asymmetry $\epsilon$")
    ax.set_ylabel("relative difference")
    ax.set_title("(d) Continuity: asymmetric code path vs symmetric path\nas the asymmetry ε → 0")
    ax.set_ylim(1e-19, 1e3)   # headroom for the annotation
    tt = extra.get("T3", {})
    if "T3_nu_int" in tt:
        fr = tt["T3_nu_int"]["frac_dev_gt_1e-6"]
        ax.text(0.03, 0.98, "values below 1e-19 drawn at 1e-19 (exact zeros)\n"
                "int, $\\nu$: fraction of entries with relative\n"
                f"  difference > 1e-6 at ε = 1e-4 / 1e-3 / 1e-2:\n  {100*fr[1]:.1f}% / {100*fr[2]:.1f}% / {100*fr[3]:.1f}%  = interference terms\n"
                "  crossing a Pauli edge (only at omega <~ 1 MeV,\n  hole windows <~ 0.5 MeV wide; set shrinks as ε → 0)",
                transform=ax.transAxes, va="top", fontsize=9)
else:
    no_data(ax)
fig.suptitle(r"B. Charge-symmetry test of the hadron tensor: $\bar\nu$ must equal $\nu$ in the mirror channel", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.97))
cap = ("Hadron-tensor unit tests: (a,b) symmetric matter, nubar tensor vs nu tensor (same mirror channel) for the new and old code; "
       "(c) asymmetric gas, nubar on the mirrored (kFn,kFp) vs nu on (kFp,kFn), grey = negative control with an un-mirrored gas (must fail); "
       "(d) continuity: relative difference between the asymmetric code path and the symmetric path as the asymmetry eps goes to 0 (median and 90th percentile). "
       "Colours: pp (nu) / nn (nubar), np (nu) / pn (nubar), pn (nu) / np (nubar).")
save_fig(fig, "B_mirror", cap, extra)
finish("B")
