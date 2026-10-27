"""E. C12 lepton-side distributions per final-state channel, flux-averaged over the flat E_nu flux (threshold to 2 GeV):
dsigma/dq0, dsigma/dE_mu, dsigma/dcos(theta_mu), linear (top) and log (bottom). One figure per beam."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = str(__import__("pathlib").Path(__file__).resolve().parents[2])  # validation/asymmetric-nuclei
COL = ["#D55E00", "#0072B2", "#009E73"]
LAB = {"nu": ["pp", "np", "pn"], "nubar": ["nn", "pn", "np"]}
SYM = {"nu": "ν", "nubar": "ν̄"}
VARS = [("q0", np.arange(0, 1.2001, 0.025), "q$_0$ = E$_\\nu$ − E$_\\mu$ [GeV]", "dσ/dq$_0$ [10$^{-38}$ cm$^2$/GeV]"),
        ("Emu", np.arange(0, 2.0001, 0.05), "E$_\\mu$ [GeV]", "dσ/dE$_\\mu$ [10$^{-38}$ cm$^2$/GeV]"),
        ("cos", np.arange(-1, 1.0001, 0.05), "cos θ$_\\mu$", "dσ/dcos θ$_\\mu$ [10$^{-38}$ cm$^2$]")]

for beam in ["nu", "nubar"]:
    d = np.load(f"{S}/data/flat_C12_{beam}_RD0.npz")
    c = d["constants"]; conv = c[3] / 1000 * 2 * np.pi / c[5]
    E = d["energy_neutrino"].astype(float); q0 = d["q0"].astype(float)
    x = {"q0": q0 / 1000, "Emu": (E - q0) / 1000, "cos": np.cos(d["theta_lepton"].astype(float))}
    w = [d[f"w{i}"] for i in range(3)]
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    for j, (key, b, xl, yl) in enumerate(VARS):
        cen = 0.5 * (b[1:] + b[:-1]); dx = np.diff(b)
        tot = np.histogram(x[key], b, weights=sum(w))[0] * conv / dx
        tote = np.sqrt(np.histogram(x[key], b, weights=sum(w) ** 2)[0]) * conv / dx
        for row in range(2):
            a = ax[row, j]
            a.step(cen, tot, where="mid", color="k", lw=1.8, label="total")
            a.fill_between(cen, tot - tote, tot + tote, step="mid", color="k", alpha=0.15, lw=0)
            for k in range(3):
                h = np.histogram(x[key], b, weights=w[k])[0] * conv / dx
                he = np.sqrt(np.histogram(x[key], b, weights=w[k] ** 2)[0]) * conv / dx
                a.step(cen, h, where="mid", color=COL[k], lw=1.5, label=LAB[beam][k])
                a.fill_between(cen, h - he, h + he, step="mid", color=COL[k], alpha=0.2, lw=0)
            a.set_xlabel(xl); a.set_ylabel(yl); a.grid(alpha=0.3)
            if row == 1:
                a.set_yscale("log"); pos = tot[tot > 0]
                a.set_ylim(max(pos.min(), pos.max() * 1e-4) * 0.5, pos.max() * 2)
            else:
                a.set_ylim(bottom=0)
        ax[0, j].set_title(yl.split(" [")[0] + " (linear)", fontsize=11); ax[1, j].set_title(yl.split(" [")[0] + " (log)", fontsize=11)
    ax[0, 0].legend(fontsize=10)
    fig.suptitle(f"E. C12 {SYM[beam]} 2p2h, Joanna's branch (FERMI off), flux-averaged over a flat flux up to 2 GeV; bands = MC statistical error", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(f"{S}/plots/figs/E_kin_C12_{beam}.png", dpi=150); plt.close(fig)
    print("E written", beam)
