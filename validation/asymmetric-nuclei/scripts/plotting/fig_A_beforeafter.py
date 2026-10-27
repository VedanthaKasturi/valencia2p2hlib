"""A. Total 2p2h cross section vs E_nu: old branch (feature/NewNuclei) vs new branch (asymmetric-nuclei), both unmodified
physics, same RNG seed and same Q, so the events are the same and only the code differs. Neutrinos only (the old branch
cannot save antineutrino events)."""
import glob, numpy as np, uproot, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = str(__import__("pathlib").Path(__file__).resolve().parents[2])  # validation/asymmetric-nuclei
R = f"{S}/runs/beforeafter"
OLD, NEW = "#8a8f98", "#2c55c7"


def load(d):
    f = uproot.open(glob.glob(f"{R}/{d}/events*.root")[0])
    t = f["treeout"]
    c = f["constants"].member("fElements")
    w = sum(t[f"e/{k}"].array(library="np") for k in ["weight_pp", "weight_np", "weight_pn"])
    return dict(E=t["e/energy_neutrino"].array(library="np"), w=w, conv=c[3] / 1000 * 2 * np.pi, nt=c[5], elo=c[2] + c[4])


edges = np.arange(150, 2001, 50.0)
cen = 0.5 * (edges[1:] + edges[:-1]) / 1000
fig, axes = plt.subplots(2, 2, figsize=(13, 7.5), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1]))
for j, (nuc, A) in enumerate([("C12", 12), ("O16", 16)]):
    o, n = load(f"old_{nuc}"), load(f"new_{nuc}")
    assert np.array_equal(o["E"], n["E"]) and o["nt"] == n["nt"], "events not aligned"
    idx = np.digitize(o["E"], edges) - 1
    m = (idx >= 0) & (idx < len(cen))
    norm = o["conv"] / (o["nt"] * np.diff(edges) / (2000 - o["elo"]))
    So = np.bincount(idx[m], o["w"][m], len(cen)) * norm
    Sn = np.bincount(idx[m], n["w"][m], len(cen)) * norm
    eo = np.sqrt(np.bincount(idx[m], o["w"][m] ** 2, len(cen))) * norm
    en = np.sqrt(np.bincount(idx[m], n["w"][m] ** 2, len(cen))) * norm
    r = Sn / So
    # same events in both samples: error of the ratio from the per-event residual w_new - r*w_old
    res = n["w"][m] - r[idx[m]] * o["w"][m]
    er = np.sqrt(np.bincount(idx[m], res ** 2, len(cen))) * norm / So
    ax = axes[0, j]
    ax.plot(cen, So, color=OLD, lw=2, label="Old code")
    ax.fill_between(cen, So - eo, So + eo, color=OLD, alpha=0.25, lw=0)
    ax.plot(cen, Sn, color=NEW, lw=2, label="Joanna's branch (FERMI off)")
    ax.fill_between(cen, Sn - en, Sn + en, color=NEW, alpha=0.25, lw=0)
    ax.set_title(f"{nuc} ν, total 2p2h (pp + np + pn)", fontsize=12)
    ax.set_ylabel("σ per nucleus [10$^{-38}$ cm$^2$]")
    ax.grid(alpha=0.3); ax.legend(fontsize=9, loc="lower right")
    tot = (n["w"][m].sum() / o["w"][m].sum())
    ax2 = axes[1, j]
    ax2.plot(cen, r, color="k", lw=1.8)
    ax2.fill_between(cen, r - er, r + er, color="k", alpha=0.2, lw=0)
    ax2.axhline(1, color="grey", lw=0.8)
    ax2.set_ylabel("Joanna's branch / old code"); ax2.set_xlabel("E$_\\nu$ [GeV]"); ax2.grid(alpha=0.3)
    ax2.set_ylim(0.95, 1.2)
    print(nuc, "flux-averaged new/old", round(tot, 4), "range of bin ratios", np.round(np.nanmin(r), 3), np.round(np.nanmax(r), 3))
fig.suptitle("A. Old code vs Joanna's branch (FERMI off): same events, same Q; bands = MC statistical error", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(f"{S}/plots/figs/A_regression.png", dpi=150)
print("written")
