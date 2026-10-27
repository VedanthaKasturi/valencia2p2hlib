"""Figure H: asymmetric isotopes: per-channel sigma ratios (H_isotopes) and area-normalised shapes (H_shapes)."""
from common import *
from matplotlib.lines import Line2D

PAIRS = [("O18", "O16"), ("Ca44", "Ca40"), ("Ar40", "Ca40")]
SHAPE_PAIRS = [("O16", "O18"), ("Ca40", "Ca44")]
need = sorted({n for p in PAIRS + SHAPE_PAIRS for n in p})
dens = load_dens()
L = pair_luminosities(dens)

# ---------------------------------------------------------------------------- reductions
D = {}
Q0E = np.arange(0, 1.2501, 0.025)
for nuc in need:
    for beam in BEAMS:
        s = load_sample(nuc, beam)
        if s is None:
            continue
        d = {"edges": s.ebins(), "Q": s.Q, "tot": [sigma_window(s, k)[0] for k in range(3)]}
        d["E"] = [s.sigma_E(s.w(k), d["edges"]) for k in range(3)]
        d["Etot"] = s.sigma_E(s.w(), d["edges"])
        d["q0"] = s.dsig(s.q0 * 1e-3, s.w(), Q0E)[0]
        d["omega"] = s.dsig(s.omega * 1e-3, s.w(), Q0E)[0]
        d["omega_min"] = float(s.omega.min()); d["q0_min"] = float(s.q0.min())
        D[(nuc, beam)] = d
        s.free("w0", "w1", "w2", "wtot", "energy_neutrino", "q0")
        del s
miss = missing_list(need)
extra = {"ratios_flux_avg": {}, "luminosity_ratios": {}, "Q": {f"{n}_{b}": d["Q"] for (n, b), d in D.items()}}

# ---------------------------------------------------------------------------- H_isotopes
fig, axs = plt.subplots(2, 3, figsize=(16, 9.5))
for i, beam in enumerate(BEAMS):
    for j, (n1, n2) in enumerate(PAIRS):
        ax = axs[i, j]
        d1, d2 = D.get((n1, beam)), D.get((n2, beam))
        ax.axhline(1, color="gray", lw=0.6)
        ax.set_xlim(0.1, 2.05)
        ax.set_title(f"{n1} / {n2}, {'$\\nu$' if beam == 'nu' else '$\\bar\\nu$'}", fontsize=11)
        ax.set_xlabel(r"$E_\nu$ [GeV]")
        if j == 0:
            ax.set_ylabel(r"$\sigma(" + n1 + r")/\sigma(" + n2 + ")$ per nucleus")
        if d1 is None or d2 is None or not L:
            no_data(ax, "samples missing: " + ", ".join(f"{n}_{beam}" for n, d in ((n1, d1), (n2, d2)) if d is None))
            continue
        edges = common_edges([d1["edges"], d2["edges"]])
        n = len(edges) - 1
        i1, i2 = np.searchsorted(d1["edges"], edges[0]), np.searchsorted(d2["edges"], edges[0])
        xc = centres(edges) * 1e-3
        series = [(SLOT_FINAL[beam][k], SLOT_COL[k], d1["E"][k][0][i1:i1 + n], d1["E"][k][1][i1:i1 + n], d2["E"][k][0][i2:i2 + n], d2["E"][k][1][i2:i2 + n], slot_pair(beam, k))
                  for k in range(3)]
        series.append(("total", TOTAL_COL, d1["Etot"][0][i1:i1 + n], d1["Etot"][1][i1:i1 + n], d2["Etot"][0][i2:i2 + n], d2["Etot"][1][i2:i2 + n], "tot"))
        for lab, col, a, ea, b, eb, pair in series:
            ok = b > 0.02 * b.max()
            r, er = ratio_err(a, ea, b, eb)
            r = np.where(ok, r, np.nan); er = np.where(ok, er, np.nan)
            ax.plot(xc, r, color=col, lw=1.8 if lab == "total" else 1.5, label=lab)
            ax.fill_between(xc, r - er, r + er, color=col, alpha=0.18, lw=0)
            Lr = L[n1][pair] / L[n2][pair]
            key = f"{n1}/{n2}_{beam}_{'tot' if lab == 'total' else lab}"
            extra["luminosity_ratios"][f"{n1}/{n2}_{pair}"] = Lr
        # flux-averaged ratio of totals per slot
        extra["ratios_flux_avg"][f"{n1}/{n2}_{beam}"] = {"slots": [d1["tot"][k] / d2["tot"][k] for k in range(3)],
                                                       "total": sum(d1["tot"]) / sum(d2["tot"]),
                                                       "L_ratio_slots": [L[n1][slot_pair(beam, k)] / L[n2][slot_pair(beam, k)] for k in range(3)],
                                                       "L_ratio_tot": L[n1]["tot"] / L[n2]["tot"]}
        ax.legend(fontsize=9, loc="best", ncol=2)
stamp_missing(fig, miss)
fig.suptitle("H. Isotope cross-section ratios per channel vs $E_\\nu$ (per nucleus); bands = MC statistical error", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.96))
save_fig(fig, "H_isotopes", "Per-channel cross-section ratios O18/O16, Ca44/Ca40, Ar40/Ca40 vs E_nu for nu (top) and nubar (bottom); dashed horizontal lines are the pair-luminosity ratios (initial pair matched to each channel) from the density profiles; flux-averaged ratios in summary.json use the common 0.15-2 GeV window.", extra)

# ---------------------------------------------------------------------------- H_shapes
fig, axs = plt.subplots(3, 4, figsize=(18, 11.5), gridspec_kw={"height_ratios": [2, 2, 1.1]})
xc = centres(Q0E)
extra2 = {}
for jj, (pa, pb) in enumerate(SHAPE_PAIRS):
    for ib, beam in enumerate(BEAMS):
        col = jj * 2 + ib
        have = [(n, D[(n, beam)]) for n in (pa, pb) if (n, beam) in D]
        for row, var in enumerate(["q0", "omega"]):
            ax = axs[row, col]
            for n, d in have:
                y = d[var]
                area = (y * np.diff(Q0E)).sum()
                ax.plot(xc, y / area, color=NUC_COL[n], label=f"{n}  (Q = {d['Q']:.2f} MeV)")
            ax.set_yscale("log")
            ax.set_ylim(1e-4, None)
            ax.set_xlabel(r"$q_0$ [GeV]" if var == "q0" else r"$\omega=q_0-Q$ [GeV]")
            ax.set_ylabel((r"area-normalised d$\sigma$/d$q_0$" if var == "q0" else r"area-normalised d$\sigma$/d$\omega$") + " [GeV$^{-1}$]")
            ax.set_title(f"{pa} vs {pb}, {'$\\nu$' if beam == 'nu' else '$\\bar\\nu$'}", fontsize=10.5)
            if have:
                ax.legend(fontsize=9, loc="upper right")
            else:
                no_data(ax)
            if var == "q0":
                for n, d in have:
                    ax.axvline(d["Q"] * 1e-3, color=NUC_COL[n], ls=":", lw=1)
        ax = axs[2, col]
        if len(have) == 2:
            (n1, d1), (n2, d2) = have
            y1 = d1["omega"] / (d1["omega"] * np.diff(Q0E)).sum(); y2 = d2["omega"] / (d2["omega"] * np.diff(Q0E)).sum()
            ok = (y1 > 1e-3 * y1.max()) & (y2 > 1e-3 * y2.max())
            with np.errstate(divide="ignore", invalid="ignore"):
                r = np.where(ok, y2 / y1, np.nan)
            ax.plot(xc, r, color=NUC_COL[n2])
            ax.axhline(1, color="k", lw=0.7)
            ax.set_ylim(0.6, 1.4)
            ax.set_ylabel(f"{n2} / {n1}\n(shape in $\\omega$)")
            extra2[f"{pa}_{pb}_{beam}"] = {"omega": xc, "ratio_shape": r, "Q": {n1: d1["Q"], n2: d2["Q"]}}
        else:
            no_data(ax)
        ax.set_xlabel(r"$\omega$ [GeV]")
stamp_missing(fig, miss)
fig.suptitle(r"H. Spectral shapes (area normalised): $q_0$ includes the Q-value shift (top), $\omega=q_0-Q$ removes the rigid shift (middle); bottom: ratio of the $\omega$ shapes", fontsize=12)
fig.tight_layout(rect=(0, 0.01, 1, 0.96))
save_fig(fig, "H_shapes", "Area-normalised dsigma/dq0 and dsigma/domega (flat flux) for O16 vs O18 and Ca40 vs Ca44, nu and nubar, with Q values in the legend; bottom row: ratio of the omega shapes (hadron-tensor change after the rigid Q shift is removed).", {"shape_ratio": extra2})
finish("H")
