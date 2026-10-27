"""H (absolute). Isotope spectra dsigma/dq0 and dsigma/domega (omega = q0 - Q), absolute (per nucleus, not normalised),
linear scale, flux-averaged over the common window E_nu = 0.15-2 GeV. Joanna's branch, Re Sigma_Delta = 0."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import NUC_COL, DATA, FIGS

SYM = {"nu": "ν", "nubar": "ν̄"}
ELO_COMMON = 150.0
bins = np.arange(0, 1.2001, 0.025)
cen = 0.5 * (bins[1:] + bins[:-1])


def spectra(nuc, beam):
    d = np.load(DATA / f"flat_{nuc}_{beam}_RD0.npz"); c = d["constants"]
    E = d["energy_neutrino"].astype(float); m = E >= ELO_COMMON
    norm = c[3] / 1000 * 2 * np.pi / (c[5] * (2000 - ELO_COMMON) / (2000 - c[2] - c[4])) / np.diff(bins)
    w = sum(d[f"w{i}"] for i in range(3))[m]
    q0 = d["q0"].astype(float)[m] / 1000; om = q0 - c[4] / 1000
    out = {}
    for key, x in [("q0", q0), ("om", om)]:
        out[key] = np.histogram(x, bins, weights=w)[0] * norm
        out[key + "_e"] = np.sqrt(np.histogram(x, bins, weights=w ** 2)[0]) * norm
    out["Q"] = c[4]
    return out


pairs = [("O16", "O18"), ("Ca40", "Ca44")]
fig, ax = plt.subplots(3, 4, figsize=(18, 11), sharex=True, gridspec_kw=dict(height_ratios=[1.2, 1.2, 0.8]))
col = 0
for a_n, b_n in pairs:
    for beam in ["nu", "nubar"]:
        A, B = spectra(a_n, beam), spectra(b_n, beam)
        for row, key, yl in [(0, "q0", "dσ/dq$_0$ per nucleus [10$^{-38}$ cm$^2$/GeV]"), (1, "om", "dσ/dω per nucleus [10$^{-38}$ cm$^2$/GeV]")]:
            x = ax[row, col]
            for n, s in [(a_n, A), (b_n, B)]:
                x.step(cen, s[key], where="mid", color=NUC_COL[n], lw=1.7, label=f"{n} (Q = {s['Q']:.2f} MeV)")
                x.fill_between(cen, s[key] - s[key + "_e"], s[key] + s[key + "_e"], step="mid", color=NUC_COL[n], alpha=0.2, lw=0)
            x.set_ylabel(yl, fontsize=9); x.grid(alpha=0.3); x.set_ylim(bottom=0); x.legend(fontsize=9)
            x.set_title(f"{a_n} vs {b_n}, {SYM[beam]}: " + ("dσ/dq$_0$" if key == "q0" else "dσ/dω,  ω = q$_0$ − Q"), fontsize=11)
        r = B["om"] / A["om"]; re_ = r * np.sqrt((B["om_e"] / B["om"]) ** 2 + (A["om_e"] / A["om"]) ** 2)
        ok = (A["om"] > 0) & (B["om"] > 0)
        ax[2, col].plot(cen[ok], r[ok], color=NUC_COL[b_n], lw=1.7)
        ax[2, col].fill_between(cen[ok], (r - re_)[ok], (r + re_)[ok], color=NUC_COL[b_n], alpha=0.2, lw=0)
        ax[2, col].axhline(1, color="grey", lw=0.8); ax[2, col].set_ylim(0.5, 1.6); ax[2, col].grid(alpha=0.3)
        ax[2, col].set_ylabel(f"dσ/dω: {b_n} / {a_n}"); ax[2, col].set_xlabel("q$_0$ or ω [GeV]")
        col += 1
fig.suptitle("H. Isotope spectra, absolute (per nucleus), linear scale; flux-averaged over E = 0.15–2 GeV; Joanna's branch (FERMI off); bands = MC statistical error", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(FIGS / "H_spectra_abs.png", dpi=150); plt.close(fig)
print("H_abs written")
