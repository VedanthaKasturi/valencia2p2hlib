"""Figure D: isospin-channel fractions vs E_nu and flux averaged."""
from common import *

fig = plt.figure(figsize=(15, 16))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.05], hspace=0.38, wspace=0.18)
axa, axb, axb2, axc = fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[2, :])
extra = {"flux_avg_fraction": {}, "dominant_fraction_E": {}, "dominant_fraction_E_100MeV": {}}
EDGES100 = np.arange(150.0, E_MAX + 1e-9, 100.0)   # 100 MeV bins starting at 150 MeV (panels b1, b2 only)
xc100 = centres(EDGES100) * 1e-3
avg_frac = {}
for nuc in NUCS:
    for beam in BEAMS:
        s = load_sample(nuc, beam)
        if s is None:
            continue
        edges = s.ebins()
        xc = centres(edges) * 1e-3
        sig = [s.sigma_E(s.w(k), edges)[0] for k in range(3)]
        tot = sum(sig)
        with np.errstate(divide="ignore", invalid="ignore"):
            f = [x / tot for x in sig]
        tots = [s.sigma_tot(k) for k in range(3)]
        avg_frac[(nuc, beam)] = [t / sum(tots) for t in tots]
        extra["flux_avg_fraction"][f"{nuc}_{beam}"] = avg_frac[(nuc, beam)]
        extra["dominant_fraction_E"][f"{nuc}_{beam}"] = {"E_GeV": xc, "f": f[0]}
        if nuc == "C12":
            for k in range(3):
                axa.plot(xc, f[k], color=SLOT_COL[k], ls=BEAM_LS[beam], marker=BEAM_MARK[beam], ms=3.2,
                         label=f"{SLOT_FINAL[beam][k]} ({BEAM_LABEL[beam]})" if True else None)
        sig100 = [s.sigma_E(s.w(k), EDGES100)[0] for k in range(3)]
        with np.errstate(divide="ignore", invalid="ignore"):
            f100 = sig100[0] / sum(sig100)
        extra["dominant_fraction_E_100MeV"][f"{nuc}_{beam}"] = {"E_GeV": xc100, "f": f100}
        (axb if beam == "nu" else axb2).plot(xc100, f100, color=NUC_COL[nuc], ls=BEAM_LS[beam], marker=BEAM_MARK[beam], ms=4, lw=1.7, label=nuc)
        s.free("w0", "w1", "w2", "wtot", "energy_neutrino")
# (a)
axa.set_xlabel(r"$E_\nu$ [GeV]"); axa.set_ylabel("fraction of 2p2h cross section")
axa.set_title(r"(a) $^{12}$C: channel fractions, $\nu$ (solid) and $\bar\nu$ (dashed); mirror channels share a colour")
if axa.get_legend_handles_labels()[0]:
    axa.legend(ncol=2, loc="center right", fontsize=9, title="final state (N1'N2')", title_fontsize=9)
else:
    no_data(axa, "C12 samples missing")
axa.axhline(0, color="k", lw=0.5)
# (b1, b2)
for ax_, beam_, lab_, fin_ in ((axb, "nu", "b1", "pp"), (axb2, "nubar", "b2", "nn")):
    ax_.set_xlabel(r"$E_\nu$ [GeV]"); ax_.set_ylabel(f"{fin_} share of the 2p2h cross section")
    ax_.set_title(f"({lab_}) {BEAM_LABEL[beam_]}: {fin_} share vs $E_\\nu$ (all nuclei, 100 MeV bins)")
    if ax_.get_legend_handles_labels()[0]:
        y0_, y1_ = ax_.get_ylim(); ax_.set_ylim(y0_, y1_ + 0.3 * (y1_ - y0_))   # headroom for the legend
        ax_.legend(ncol=4, fontsize=9, loc="upper center")
    else:
        no_data(ax_)
# (c) stacked bars
if avg_frac:
    x0 = 0
    xt, xl = [], []
    gap = 0.7
    for beam in BEAMS:
        for i, nuc in enumerate(NUCS):
            x = x0 + i
            if (nuc, beam) not in avg_frac:
                axc.text(x, 0.5, "n/a", ha="center", va="center", rotation=90, color="#a00000", fontsize=9)
            else:
                f = avg_frac[(nuc, beam)]
                bot = 0.0
                for k in range(3):
                    axc.bar(x, f[k], 0.75, bottom=bot, color=SLOT_COL[k], edgecolor="white", lw=0.5,
                            label=mirror_label(k) if (x0 == 0 and i == 0) else None)
                    if abs(f[k]) > 0.04:
                        axc.text(x, bot + f[k] / 2, f"{100*f[k]:.0f}%", ha="center", va="center", fontsize=9, color="white")
                    bot += f[k]
            xt.append(x); xl.append(nuc)
        x0 += len(NUCS) + gap
    axc.set_xticks(xt); axc.set_xticklabels(xl, rotation=0)
    ymax = axc.get_ylim()[1]
    axc.text(3, 1.04, r"neutrino  (bottom to top: pp, np, pn)", ha="center", fontsize=10.5)
    axc.text(len(NUCS) + gap + 3, 1.04, r"antineutrino  (bottom to top: nn, pn, np)", ha="center", fontsize=10.5)
    axc.set_ylim(0, 1.12)
    axc.legend(loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.28), fontsize=9)
    axc.set_ylabel("flux-averaged fraction of $\\sigma$")
    axc.set_title("(c) flux-averaged channel fractions (flat flux 0-2 GeV) per nucleus and beam (stacked; percentages inside)")
    axc.grid(axis="x", visible=False)
else:
    no_data(axc)
stamp_missing(fig, missing_list())
fig.suptitle("D. Isospin-channel fractions of the 2p2h cross section", fontsize=12.5, y=0.985)
cap = ("Isospin-channel fractions: (a) 12C vs E_nu for nu (pp,np,pn; solid) and nubar (nn,pn,np; dashed), mirror channels in the same colour; "
       "(b1, b2) pp share for nu (left) and nn share for nubar (right) vs E_nu for all nuclei, 100 MeV bins starting at 150 MeV; "
       "(c) flux-averaged channel fractions per nucleus and beam.")
fig.subplots_adjust(left=0.06, right=0.985, top=0.94, bottom=0.07)
save_fig(fig, "D_fractions", cap, extra)
finish("D")
