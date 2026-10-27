"""L. Per-nucleon 2p2h cross section vs E for C12, Ar40 and Ca40, ratio to C12. Joanna's branch, Re Sigma_Delta = 0.
L1: each nucleus at its own Q.  L2: all three at the Ar40 Q (S_2p for nu, S_2n for nubar; Q2P2H_OVERRIDE runs qAr_*),
so that Ar40 vs Ca40 isolates the N/Z asymmetry and the comparison with C12 isolates the density profile."""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import A_OF, NUC_COL, DATA, FIGS

SYM = {"nu": "ν", "nubar": "ν̄"}
NUCS3 = ["C12", "Ar40", "Ca40"]
ELO_COMMON = 150.0
VERSIONS = {  # figure -> {nucleus: data tag prefix}
    "L1_C12_Ar40_Ca40_ownQ": {"C12": "flat_C12", "Ar40": "flat_Ar40", "Ca40": "flat_Ca40"},
    "L2_C12_Ar40_Ca40_commonQ": {"C12": "qAr_C12", "Ar40": "flat_Ar40", "Ca40": "qAr_Ca40"},
}
TITLE = {"L1_C12_Ar40_Ca40_ownQ": "each nucleus at its own Q",
         "L2_C12_Ar40_Ca40_commonQ": "all at the Ar40 Q"}


def load(prefix, beam):
    d = np.load(DATA / f"{prefix}_{beam}_RD0.npz")
    c = d["constants"]
    return dict(E=d["energy_neutrino"].astype(float), w=np.vstack([d[f"w{i}"] for i in range(3)]).sum(0),
                conv=c[3] / 1000 * 2 * np.pi, nt=c[5], elo=c[2] + c[4], Q=c[4])


def sigma_E(s, edges):
    idx = np.digitize(s["E"], edges) - 1; m = (idx >= 0) & (idx < len(edges) - 1)
    norm = s["conv"] / (s["nt"] * np.diff(edges) / (2000 - s["elo"]))
    return np.bincount(idx[m], s["w"][m], len(edges) - 1) * norm, np.sqrt(np.bincount(idx[m], s["w"][m] ** 2, len(edges) - 1)) * norm


def flux_avg(s):
    """flux-averaged sigma over the common window 0.15-2 GeV with its MC error"""
    m = s["E"] >= ELO_COMMON
    norm = s["conv"] / (s["nt"] * (2000 - ELO_COMMON) / (2000 - s["elo"]))
    return s["w"][m].sum() * norm, np.sqrt((s["w"][m] ** 2).sum()) * norm


edges = np.append(np.arange(150, 1951, 100.0), 2000.0)
cen = 0.5 * (edges[1:] + edges[:-1]) / 1000
summary = {}
for name, tags in VERSIONS.items():
    fig, ax = plt.subplots(2, 2, figsize=(13, 8.5), sharex=True, gridspec_kw=dict(height_ratios=[1.7, 1]))
    for j, b in enumerate(["nu", "nubar"]):
        S = {n: load(tags[n], b) for n in NUCS3}
        ref, refe = sigma_E(S["C12"], edges); ref, refe = ref / 12, refe / 12
        fref, frefe = flux_avg(S["C12"]); fref, frefe = fref / 12, frefe / 12
        for n in NUCS3:
            v, e = sigma_E(S[n], edges); v, e = v / A_OF[n], e / A_OF[n]
            f, fe = flux_avg(S[n]); f, fe = f / A_OF[n], fe / A_OF[n]
            fr = f / fref; fre = 0.0 if n == "C12" else fr * np.hypot(fe / f, frefe / fref)
            summary.setdefault(name, {})[f"{n}_{b}"] = dict(Q=float(S[n]["Q"]), sigmaA_fluxavg=f, err=fe,
                                                         ratio_to_C12=fr, ratio_err=fre)
            lab = f"{n} (Q = {S[n]['Q']:.1f} MeV)" + ("" if n == "C12" else f",  ⟨σ/A⟩/C12 = {fr:.3f} ± {fre:.3f}")
            ax[0, j].plot(cen, v, color=NUC_COL[n], lw=1.8, label=lab)
            ax[0, j].fill_between(cen, v - e, v + e, color=NUC_COL[n], alpha=0.18, lw=0)
            if n != "C12":
                r = v / ref; re_ = r * np.sqrt((e / v) ** 2 + (refe / ref) ** 2)
                ax[1, j].plot(cen, r, color=NUC_COL[n], lw=1.8)
                ax[1, j].fill_between(cen, r - re_, r + re_, color=NUC_COL[n], alpha=0.18, lw=0)
        ax[0, j].set_title(f"{SYM[b]}: 2p2h cross section per nucleon", fontsize=12)
        ax[0, j].set_ylabel("σ / A  [10$^{-38}$ cm$^2$ per nucleon]"); ax[0, j].grid(alpha=0.3); ax[0, j].legend(fontsize=9)
        ax[1, j].axhline(1, color="grey", lw=0.8); ax[1, j].set_ylabel("(σ/A) / (σ/A)$_{C12}$")
        ax[1, j].set_xlabel(f"E$_{{{SYM[b]}}}$ [GeV]"); ax[1, j].grid(alpha=0.3)
    fig.suptitle(f"{name.split('_')[0]}. Per-nucleon 2p2h cross section, C12 / Ar40 / Ca40, {TITLE[name]}; "
                 "100 MeV bins, MC statistical bands; ⟨⟩ = flux average over 0.15–2 GeV", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIGS / f"{name}.png", dpi=150); plt.close(fig)
    print("wrote", FIGS / f"{name}.png")

(FIGS / "L_C12_Ar40_Ca40.json").write_text(json.dumps(summary, indent=1))
for name, d in summary.items():
    print(name)
    for k, v in d.items():
        print(f"  {k:12s} Q={v['Q']:7.3f}  <σ/A>={v['sigmaA_fluxavg']:.5f}  /C12 = {v['ratio_to_C12']:.3f} ± {v['ratio_err']:.3f}")
