"""C. C12 antineutrino 2p2h from the generated events (flat E_nu flux), 100 MeV E_nu bins.
C1: old code vs Joanna's branch, same Q = 27.41 MeV: total and the nn / pn / np contributions.
C2: Joanna's branch as delivered (Re Sigma_Delta = 0): total and channels.
C3: Joanna's branch vs PRC 111 Fig. 8, with the Re Sigma_Delta band."""
import glob, numpy as np, uproot, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = str(__import__("pathlib").Path(__file__).resolve().parents[2])  # validation/asymmetric-nuclei
FIG = f"{S}/plots/figs"
CH = ["nn", "pn", "np"]
COL = ["#D55E00", "#0072B2", "#009E73"]
edges = np.append(np.arange(150, 1951, 100.0), 2000.0)  # centres 0.2 ... 1.9 GeV (one on 1 GeV) plus a last 1.95-2.00 GeV bin
cen = 0.5 * (edges[1:] + edges[:-1]) / 1000


def from_root(d):
    f = uproot.open(glob.glob(f"{S}/runs/oldnubar/{d}/events*.root")[0])
    t = f["treeout"]; c = f["constants"].member("fElements")
    w = np.vstack([t[f"e/{k}"].array(library="np") for k in ["weight_anti_nn", "weight_anti_pn", "weight_anti_np"]])
    return dict(E=t["e/energy_neutrino"].array(library="np"), w=w, c=c)


def from_npz(tag):
    d = np.load(f"{S}/data/{tag}.npz")
    return dict(E=d["energy_neutrino"].astype(float), w=np.vstack([d[f"w{i}"] for i in range(3)]), c=d["constants"])


def sigma(s):
    """sigma(E) per channel [3, nbins] and MC error, 1e-38 cm^2 per nucleus."""
    c = s["c"]; conv = c[3] / 1000 * 2 * np.pi; elo = c[2] + c[4]
    idx = np.digitize(s["E"], edges) - 1; m = (idx >= 0) & (idx < len(cen))
    norm = conv / (c[5] * np.diff(edges) / (2000 - elo))
    v = np.vstack([np.bincount(idx[m], s["w"][k][m], len(cen)) * norm for k in range(3)])
    e = np.vstack([np.sqrt(np.bincount(idx[m], s["w"][k][m] ** 2, len(cen))) * norm for k in range(3)])
    et = np.sqrt(np.bincount(idx[m], s["w"].sum(0)[m] ** 2, len(cen))) * norm
    return v, e, et


def stacked(ax, v, et, title):
    base = np.zeros_like(cen)
    for k in range(3):
        ax.fill_between(cen, base, base + v[k], step="mid", color=COL[k], alpha=0.75, lw=0, label=CH[k])
        base = base + v[k]
    ax.step(cen, base, where="mid", color="k", lw=1.6, label="total")
    ax.fill_between(cen, base - et, base + et, step="mid", color="k", alpha=0.15, lw=0)
    ax.set_title(title, fontsize=12); ax.set_ylabel("σ per nucleus [10$^{-38}$ cm$^2$]"); ax.grid(alpha=0.3)


def fractions(ax, v, e):
    tot = v.sum(0)
    for k in range(3):
        f = v[k] / tot
        ax.plot(cen, 100 * f, color=COL[k], lw=1.8, marker="o", ms=3, label=CH[k])
        ax.fill_between(cen, 100 * (f - e[k] / tot), 100 * (f + e[k] / tot), color=COL[k], alpha=0.2, lw=0)
    ax.set_ylabel("share of total [%]"); ax.set_xlabel("E$_{\\bar\\nu}$ [GeV]"); ax.set_ylim(0, 100); ax.grid(alpha=0.3)


# ---------- C1: old code vs Joanna's branch, same Q
old, new = from_root("old_C12_nubar"), from_root("new_C12_nubar_Q27p41")
vo, eo, eto = sigma(old); vn, en, etn = sigma(new)
fig, ax = plt.subplots(2, 2, figsize=(13, 8.5), sharex=True, sharey="row", gridspec_kw=dict(height_ratios=[1.6, 1]))
stacked(ax[0, 0], vo, eto, "Old code"); stacked(ax[0, 1], vn, etn, "Joanna's branch (FERMI off)")
fractions(ax[1, 0], vo, eo); fractions(ax[1, 1], vn, en)
ax[0, 0].legend(fontsize=10, loc="upper left"); ax[1, 1].legend(fontsize=10, loc="center right")
fig.suptitle("C1. C12 ν̄ 2p2h: total and its nn / pn / np contributions (same Q = 27.41 MeV in both)", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(f"{FIG}/C1_nubar_old_vs_new.png", dpi=150); plt.close(fig)
fo, fn = vo.sum(1) / vo.sum(), vn.sum(1) / vn.sum()
print("C1 shares old:", np.round(100 * fo, 1), " new:", np.round(100 * fn, 1), " total new/old:", round(vn.sum() / vo.sum(), 3))

# ---------- C2: Joanna's branch as delivered, Re Sigma_Delta = 0
s0 = from_npz("flat_C12_nubar_RD0"); v0, e0, et0 = sigma(s0)
fig, ax = plt.subplots(2, 1, figsize=(8.5, 8.5), sharex=True, gridspec_kw=dict(height_ratios=[1.6, 1]))
stacked(ax[0], v0, et0, f"Joanna's branch (FERMI off), Q = {s0['c'][4]:.2f} MeV, Re Σ$_\\Delta$ = 0")
fractions(ax[1], v0, e0); ax[0].legend(fontsize=10, loc="upper left")
fig.suptitle("C2. C12 ν̄ 2p2h: total and its nn / pn / np contributions", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(f"{FIG}/C2_nubar_channels.png", dpi=150); plt.close(fig)

# ---------- C3: comparison with PRC 111 Fig. 8, with the Re Sigma_Delta band
s1, s2 = from_npz("flat_C12_nubar_RD1"), from_npz("flat_C12_nubar_RD2")
t0, t1, t2 = v0.sum(0), sigma(s1)[0].sum(0), sigma(s2)[0].sum(0)
lo, hi = np.minimum(t1, t2), np.maximum(t1, t2)
paper = np.array([(0.25, 0.0), (0.375, 0.09), (0.5, 0.2), (0.625, 0.33), (0.75, 0.45), (0.875, 0.55), (1.0, 0.63),
                  (1.25, 0.78), (1.5, 0.9), (1.75, 0.99), (2.0, 1.13)])
perr = np.maximum(0.05 * paper[:, 1], 0.05)
fig, ax = plt.subplots(2, 1, figsize=(9.5, 8.5), sharex=True, gridspec_kw=dict(height_ratios=[1.8, 1]))
a = ax[0]
a.fill_between(cen, lo, hi, step="mid", color="#2c55c7", alpha=0.2, lw=0, label="Joanna's branch, Re Σ$_\\Delta$ = −50 … +30 ρ/ρ$_0$ MeV")
a.step(cen, t0, where="mid", color="#2c55c7", lw=2, label="Joanna's branch, Re Σ$_\\Delta$ = 0 (2p2h only)")
a.fill_between(cen, t0 - et0, t0 + et0, step="mid", color="#2c55c7", alpha=0.35, lw=0)
a.errorbar(paper[:, 0], paper[:, 1], yerr=perr, fmt="o", color="#D55E00", ms=5, capsize=2, label="PRC 111 Fig. 8, central curve (read by eye)")
a.errorbar([1.0], [0.625], yerr=[[0.085], [0.085]], fmt="none", ecolor="#D55E00", elinewidth=3, capsize=7, label="PRC 111 text at 1 GeV: band [0.54, 0.71]")
a.set_ylabel("σ per nucleus [10$^{-38}$ cm$^2$]"); a.grid(alpha=0.3); a.legend(fontsize=9, loc="upper left")
a.set_title("C12 ν̄, FERMI off. The paper curves include 3p3h; it is not added here", fontsize=11)
r = np.interp(cen, paper[:, 0], paper[:, 1]); ok = (cen >= 0.375) & (cen <= 2.0)  # below 0.375 GeV the paper value is ~0
ax[1].plot(cen[ok], (t0 / r)[ok], color="#2c55c7", lw=2, marker="o", ms=3, label="Re Σ$_\\Delta$ = 0")
ax[1].fill_between(cen[ok], (lo / r)[ok], (hi / r)[ok], color="#2c55c7", alpha=0.2, lw=0)
ax[1].axhline(1, color="grey", lw=0.8); ax[1].set_ylim(0.5, 1.3)
ax[1].set_ylabel("this 2p2h / paper (2p2h + 3p3h)")
ax[1].text(0.01, 0.04, "ratio from 0.4 GeV: below that the paper value is ~0", transform=ax[1].transAxes, fontsize=9); ax[1].set_xlabel("E$_{\\bar\\nu}$ [GeV]"); ax[1].grid(alpha=0.3)
fig.suptitle("C3. C12 ν̄ 2p2h, Joanna's branch vs Sobczyk & Nieves PRC 111, 025502, Fig. 8", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(f"{FIG}/C3_nubar_vs_paper.png", dpi=150); plt.close(fig)
i1 = np.argmin(np.abs(cen - 1.0))
print(f"C3 at {cen[i1]:.2f} GeV: RD0 {t0[i1]:.3f}  band [{lo[i1]:.3f}, {hi[i1]:.3f}]")
