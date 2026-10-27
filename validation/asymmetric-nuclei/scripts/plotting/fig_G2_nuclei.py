"""G (revised). Per-nucleon 2p2h cross section vs E for all nuclei (with and without MC bands), and scaling of the
flux-averaged cross section with the neutron number N and the proton number Z. Joanna's branch, Re Sigma_Delta = 0."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import NUCS, A_OF, Z_OF, N_OF, NUC_COL, SLOT_COL, DATA, FIGS

LAB = {"nu": ["pp", "np", "pn"], "nubar": ["nn", "pn", "np"]}
SYM = {"nu": "ν", "nubar": "ν̄"}
ELO_COMMON = 150.0


def load(nuc, beam):
    d = np.load(DATA / f"flat_{nuc}_{beam}_RD0.npz")
    c = d["constants"]
    return dict(E=d["energy_neutrino"].astype(float), w=np.vstack([d[f"w{i}"] for i in range(3)]),
                conv=c[3] / 1000 * 2 * np.pi, nt=c[5], elo=c[2] + c[4], Q=c[4])


def sigma_E(s, edges):
    idx = np.digitize(s["E"], edges) - 1; m = (idx >= 0) & (idx < len(edges) - 1)
    norm = s["conv"] / (s["nt"] * np.diff(edges) / (2000 - s["elo"]))
    w = s["w"].sum(0)
    return np.bincount(idx[m], w[m], len(edges) - 1) * norm, np.sqrt(np.bincount(idx[m], w[m] ** 2, len(edges) - 1)) * norm


def flux_avg(s):
    """flux-averaged sigma over the common window 0.15-2 GeV, per channel and total, with MC errors"""
    m = s["E"] >= ELO_COMMON
    norm = s["conv"] / (s["nt"] * (2000 - ELO_COMMON) / (2000 - s["elo"]))
    ch = s["w"][:, m].sum(1) * norm; che = np.sqrt((s["w"][:, m] ** 2).sum(1)) * norm
    t = s["w"][:, m].sum(0)
    return ch, che, t.sum() * norm, np.sqrt((t ** 2).sum()) * norm


S = {(n, b): load(n, b) for n in NUCS for b in ["nu", "nubar"]}

# ---------- per-nucleon sigma(E), without and with bands
edges = np.append(np.arange(150, 1951, 100.0), 2000.0)
cen = 0.5 * (edges[1:] + edges[:-1]) / 1000
for bands in [False, True]:
    fig, ax = plt.subplots(2, 2, figsize=(14, 9), sharex=True, gridspec_kw=dict(height_ratios=[1.7, 1]))
    for j, b in enumerate(["nu", "nubar"]):
        ref, refe = sigma_E(S[("C12", b)], edges); ref, refe = ref / 12, refe / 12
        for n in NUCS:
            v, e = sigma_E(S[(n, b)], edges); v, e = v / A_OF[n], e / A_OF[n]
            ax[0, j].plot(cen, v, color=NUC_COL[n], lw=1.8, label=f"{n} (Q = {S[(n, b)]['Q']:.1f} MeV)")
            r = v / ref
            ax[1, j].plot(cen, r, color=NUC_COL[n], lw=1.8)
            if bands:
                ax[0, j].fill_between(cen, v - e, v + e, color=NUC_COL[n], alpha=0.18, lw=0)
                re_ = r * np.sqrt((e / v) ** 2 + (refe / ref) ** 2)
                ax[1, j].fill_between(cen, r - re_, r + re_, color=NUC_COL[n], alpha=0.18, lw=0)
        ax[0, j].set_title(f"{SYM[b]}: 2p2h cross section per nucleon", fontsize=12)
        ax[0, j].set_ylabel("σ / A  [10$^{-38}$ cm$^2$ per nucleon]"); ax[0, j].grid(alpha=0.3); ax[0, j].legend(fontsize=9)
        ax[1, j].axhline(1, color="grey", lw=0.8); ax[1, j].set_ylabel("(σ/A) / (σ/A)$_{C12}$")
        ax[1, j].set_xlabel(f"E$_{{{'ν' if b == 'nu' else 'ν̄'}}}$ [GeV]"); ax[1, j].grid(alpha=0.3)
    tag = "with MC statistical bands" if bands else "no error bands"
    fig.suptitle(f"G1. Per-nucleon 2p2h cross section, Joanna's branch (FERMI off), 100 MeV bins, {tag}", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIGS / f"G1_sigmaA_{'bands' if bands else 'nobands'}.png", dpi=150); plt.close(fig)

# ---------- scaling with N and Z
fig, ax = plt.subplots(2, 2, figsize=(14, 10.5))
out = {}
for j, b in enumerate(["nu", "nubar"]):
    vals = {n: flux_avg(S[(n, b)]) for n in NUCS}
    series = [("total", "k", lambda v: (v[2], v[3])),
              (LAB[b][0], SLOT_COL[0], lambda v: (v[0][0], v[1][0])),
              (f"{LAB[b][1]} + {LAB[b][2]}", SLOT_COL[1], lambda v: (v[0][1] + v[0][2], np.hypot(v[1][1], v[1][2])))]
    for i, (xname, XOF) in enumerate([("N (neutrons)", N_OF), ("Z (protons)", Z_OF)]):
        a = ax[i, j]
        x = np.array([XOF[n] for n in NUCS], float)
        for lab, col, f in series:
            y = np.array([f(vals[n])[0] for n in NUCS]); ye = np.array([f(vals[n])[1] for n in NUCS])
            slope, icpt = np.polyfit(np.log(x), np.log(y), 1)
            xx = np.linspace(x.min() * 0.9, x.max() * 1.1, 50)
            a.errorbar(x, y, yerr=ye, fmt="o", color=col, ms=6, capsize=2, label=f"{lab}: σ ∝ {xname[0]}$^{{{slope:.2f}}}$")
            a.plot(xx, np.exp(icpt) * xx ** slope, color=col, lw=0.9, ls="--", alpha=0.7)
            out[f"{b}_{lab}_vs_{xname[0]}"] = round(float(slope), 3)
        for n in NUCS:
            y = vals[n][2]
            off = {"C12": (6, -2), "O16": (-28, 4), "O18": (6, -12), "Ar40": (-30, 6), "Ca40": (6, -13), "Ca44": (6, 6), "Pb208": (-38, 6)}[n]
            a.annotate(n, (XOF[n], y), textcoords="offset points", xytext=off, fontsize=8.5, color=NUC_COL[n])
        a.set_xscale("log"); a.set_yscale("log")
        a.set_xticks([6, 8, 10, 20, 40, 80, 126]); a.set_xticklabels(["6", "8", "10", "20", "40", "80", "126"])
        a.set_xlabel(xname); a.set_ylabel("flux-averaged σ per nucleus [10$^{-38}$ cm$^2$]")
        a.set_title(f"{SYM[b]}: σ vs {xname}", fontsize=12); a.grid(alpha=0.3, which="both"); a.legend(fontsize=9, loc="upper left")
fig.suptitle("G2. Flux-averaged (0.15–2 GeV) 2p2h cross section vs neutron and proton number; dashed: power-law fit over the 7 nuclei", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(FIGS / "G2_scaling_ZN.png", dpi=150); plt.close(fig)
print("exponents:", out)

# ---------- scaling with A
fig, ax = plt.subplots(2, 2, figsize=(14, 10.5))
outA = {}
x = np.array([A_OF[n] for n in NUCS], float)
nz = np.array([Z_OF[n] == N_OF[n] for n in NUCS])
for j, b in enumerate(["nu", "nubar"]):
    vals = {n: flux_avg(S[(n, b)]) for n in NUCS}
    series = [("total", "k", lambda v: (v[2], v[3])),
              (LAB[b][0], SLOT_COL[0], lambda v: (v[0][0], v[1][0])),
              (f"{LAB[b][1]} + {LAB[b][2]}", SLOT_COL[1], lambda v: (v[0][1] + v[0][2], np.hypot(v[1][1], v[1][2])))]
    a = ax[0, j]
    for lab, col, f in series:
        y = np.array([f(vals[n])[0] for n in NUCS]); ye = np.array([f(vals[n])[1] for n in NUCS])
        sl, ic = np.polyfit(np.log(x), np.log(y), 1)
        xx = np.linspace(10, 230, 50)
        a.errorbar(x, y, yerr=ye, fmt="o", color=col, ms=6, capsize=2, label=f"{lab}: σ ∝ A$^{{{sl:.2f}}}$")
        a.plot(xx, np.exp(ic) * xx ** sl, color=col, lw=0.9, ls="--", alpha=0.7)
        outA[f"{b}_{lab}"] = round(float(sl), 3)
    y0 = vals["C12"][2]
    a.plot(xx, y0 * xx / 12, color="grey", lw=1.2, ls=":", label="σ ∝ A (through C12)")
    a.set_xscale("log"); a.set_yscale("log"); a.set_xticks([12, 16, 20, 40, 100, 208]); a.set_xticklabels(["12", "16", "20", "40", "100", "208"])
    a.set_xlabel("A = Z + N"); a.set_ylabel("flux-averaged σ per nucleus [10$^{-38}$ cm$^2$]")
    a.set_title(f"{SYM[b]}: σ vs A", fontsize=12); a.grid(alpha=0.3, which="both"); a.legend(fontsize=9, loc="upper left")
    # per nucleon vs A, N = Z filled, N > Z open
    a = ax[1, j]
    y = np.array([vals[n][2] / A_OF[n] for n in NUCS]); ye = np.array([vals[n][3] / A_OF[n] for n in NUCS])
    sl, ic = np.polyfit(np.log(x), np.log(y * x), 1)
    a.plot(xx, np.exp(ic) * xx ** (sl - 1), color="k", lw=0.9, ls="--", alpha=0.7, label=f"fit: σ/A ∝ A$^{{{sl - 1:.2f}}}$")
    for i, n in enumerate(NUCS):
        a.errorbar(x[i], y[i], yerr=ye[i], fmt="o", ms=8, capsize=2, color=NUC_COL[n],
                   mfc=NUC_COL[n] if nz[i] else "white", mew=1.8)
        off = {"C12": (8, -3), "O16": (-30, 6), "O18": (8, -12), "Ar40": (-34, 6), "Ca40": (8, -12), "Ca44": (8, 4), "Pb208": (-40, 8)}[n]
        a.annotate(n, (x[i], y[i]), textcoords="offset points", xytext=off, fontsize=9, color=NUC_COL[n])
    a.errorbar([], [], fmt="o", color="grey", mfc="grey", label="N = Z"); a.errorbar([], [], fmt="o", color="grey", mfc="white", label="N > Z")
    a.set_xscale("log"); a.set_xticks([12, 16, 20, 40, 100, 208]); a.set_xticklabels(["12", "16", "20", "40", "100", "208"])
    a.set_xlabel("A = Z + N"); a.set_ylabel("flux-averaged σ / A [10$^{-38}$ cm$^2$ per nucleon]")
    a.set_title(f"{SYM[b]}: σ per nucleon vs A", fontsize=12); a.grid(alpha=0.3, which="both"); a.legend(fontsize=9, loc="upper left")
fig.suptitle("G3. Flux-averaged (0.15–2 GeV) 2p2h cross section vs mass number A; dashed: power-law fit over the 7 nuclei", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(FIGS / "G3_scaling_A.png", dpi=150); plt.close(fig)
print("A exponents:", outA)
