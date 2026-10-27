"""Figure I: radial dependence: density profiles / kF ratio, dsigma/dr, channel fractions vs r."""
from common import *
from matplotlib.lines import Line2D

dens = load_dens()
RP_COL, RN_COL = "#CC3311", "#0077BB"
DR = 0.2
FR_NUC = ["Ar40", "Ca40", "Pb208"]
RD = {}
for nuc in NUCS:
    for beam in BEAMS:
        s = load_sample(nuc, beam)
        if s is None:
            continue
        r = s.a("r")
        rmax = np.ceil(r.max() / DR) * DR
        edges = np.arange(0, rmax + 1e-9, DR)
        b = Binner(r, edges)
        d = {"edges": edges}
        d["slot"] = [b.sums(s.w(k))[0] for k in range(3)]
        d["tot"] = sum(d["slot"])
        d["n"] = b.count()
        d["fac"] = s.fac
        RD[(nuc, beam)] = d
        s.free("w0", "w1", "w2", "wtot")
        del s
miss = missing_list()

fig = plt.figure(figsize=(19, 22))
gs = fig.add_gridspec(5, 4, height_ratios=[1, 1, 1.05, 1, 1], hspace=0.42, wspace=0.55)
extra = {}
# (a) densities --------------------------------------------------------------------------------------
if dens is not None:
    for i, nuc in enumerate(NUCS):
        ax = fig.add_subplot(gs[i // 4, i % 4])
        d = dens[dens.nucleus == nuc]
        ax.plot(d.r_fm, d.rho_p_fm3, color=RP_COL, lw=2.2, label=r"$\rho_p$")
        ax.plot(d.r_fm, d.rho_n_fm3, color=RN_COL, ls="--", label=r"$\rho_n$")
        ax.set_xlabel("r [fm]"); ax.set_ylabel(r"density [fm$^{-3}$]")
        ax.set_title(f"{nuc}  (Z={Z_OF[nuc]}, N={N_OF[nuc]})", fontsize=11, color=NUC_COL[nuc] if nuc != "O16" else "#2a8fb8")
        ax2 = ax.twinx()
        with np.errstate(divide="ignore", invalid="ignore"):
            kr = d.kF_n_MeV / d.kF_p_MeV
        ok = np.isfinite(kr) & (d.kF_p_MeV > 0)
        ax2.plot(d.r_fm[ok], kr[ok], color="k", ls=":", lw=1.5, label=r"$k_{F,n}/k_{F,p}$")
        ax2.set_ylim(0.9, max(1.5, float(kr[ok].max()) + 0.02)); ax2.set_ylabel(r"$k_{F,n}/k_{F,p}$"); ax2.grid(False); ax2.spines["right"].set_visible(True)
        extra[nuc] = {"kFn_over_kFp_center": float(kr[ok].iloc[0]), "kFn_over_kFp_max": float(kr[ok].max())}
        if i == 0:
            h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
            ax.legend(h1 + h2, l1 + l2, fontsize=9, loc="upper right")
    axl = fig.add_subplot(gs[1, 3]); axl.axis("off")
    axl.text(0, 0.9, "(a) density profiles from densities.csv\n(proton, neutron; kF ratio on the right axis)", fontsize=10.5, va="top")
else:
    ax = fig.add_subplot(gs[0:2, :]); no_data(ax, "densities.csv missing")
# (b) dsigma/dr normalised ------------------------------------------------------------------------------------
axb = fig.add_subplot(gs[2, :2]); axb2 = fig.add_subplot(gs[2, 2:])
for ax, beam in zip((axb, axb2), BEAMS):
    drawn = False
    for nuc in NUCS:
        d = RD.get((nuc, beam))
        if d is None:
            continue
        xc = centres(d["edges"])
        y = d["tot"] / (d["tot"].sum() * DR)
        ax.plot(xc, y, color=NUC_COL[nuc], label=nuc)
        drawn = True
        extra.setdefault(nuc, {})[f"mean_r_{beam}"] = float((xc * d["tot"]).sum() / d["tot"].sum())
    ax.set_xlabel("r [fm]"); ax.set_ylabel(r"normalised d$\sigma$/d$r$  [fm$^{-1}$]")
    ax.set_title(f"(b) {'neutrino' if beam == 'nu' else 'antineutrino'}: radial distribution of the 2p2h cross section (area = 1)", fontsize=10.5)
    if drawn:
        ax.legend(fontsize=9, ncol=2)
    else:
        no_data(ax)
# (c) slot fractions vs r --------------------------------------------------------------------------------------------
for row, beam in enumerate(BEAMS):
    for col, nuc in enumerate(FR_NUC):
        ax = fig.add_subplot(gs[3 + row, col])
        d = RD.get((nuc, beam))
        ax.set_title(f"(c) {nuc}, {'$\\nu$' if beam == 'nu' else '$\\bar\\nu$'}", fontsize=10.5)
        ax.set_xlabel("r [fm]"); ax.set_ylabel("channel fraction of $\\sigma$ in the r-bin")
        if d is None:
            no_data(ax, f"{nuc}_{beam} missing")
            continue
        # coarser bins (0.6 fm) to tame noise
        k = 3
        nb = (len(d["tot"]) // k) * k
        xc = d["edges"][:nb:k] + k * DR / 2
        tot = d["tot"][:nb].reshape(-1, k).sum(1)
        cnt = d["n"][:nb].reshape(-1, k).sum(1)
        for sl in range(3):
            y = d["slot"][sl][:nb].reshape(-1, k).sum(1) / tot
            ax.plot(xc, np.where(cnt > 500, y, np.nan), color=SLOT_COL[sl], label=SLOT_FINAL[beam][sl])
        ax.legend(fontsize=9, loc="center right")
    # right-most cell in the row: text
axt = fig.add_subplot(gs[3:5, 3]); axt.axis("off")
axt.text(0, 0.95, "(c) channel fractions vs radius\n(0.6 fm bins; bins with <500 events omitted).\n\nFinal-state channels:\npp, np, pn ($\\nu$);\nnn, pn, np ($\\bar\\nu$).", fontsize=10.5, va="top")
stamp_missing(fig, miss)
fig.suptitle("I. Radial dependence: densities and local Fermi momenta (a), radial distribution of the cross section (b), channel fractions vs r (c)", fontsize=13, y=0.992)
fig.subplots_adjust(left=0.05, right=0.975, top=0.965, bottom=0.04)
save_fig(fig, "I_radial", "Radial structure: (a) proton/neutron density profiles and kF_n/kF_p(r) per nucleus; (b) area-normalised dsigma/dr for nu and nubar, all nuclei; (c) channel fractions vs r for Ar40, Ca40 and Pb208.", {"nuclei": extra})
finish("I")
