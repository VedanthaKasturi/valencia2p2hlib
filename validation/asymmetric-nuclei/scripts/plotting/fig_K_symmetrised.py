"""K. Effect of the per-species Fermi seas alone: asymmetric treatment (branch) vs symmetrised local Fermi gas
(rho_p = rho_n = rho/2 at every r, i.e. the old single-kF treatment), same nucleus, same Q, same everything else.
Inputs: S/data/flat_<NUC>_<beam>_RD0.npz (branch) and S/data/sym_<NUC>_<beam>_RD0.npz (symmetrised control)."""
import os, json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = str(__import__("pathlib").Path(__file__).resolve().parents[2])  # validation/asymmetric-nuclei
FIG = f"{S}/plots/figs"
NUCS = ["O18", "Ar40", "Ca40", "Ca44", "Pb208"]
NZ = {"O18": (8, 10), "Ar40": (18, 22), "Ca40": (20, 20), "Ca44": (20, 24), "Pb208": (82, 126)}
SLOTCOL = ["#D55E00", "#0072B2", "#009E73"]
LAB = {"nu": ["pp", "np", "pn"], "nubar": ["nn", "pn", "np"]}


def load(tag):
    p = f"{S}/data/{tag}.npz"
    if not os.path.exists(p):
        return None
    d = np.load(p)
    c = d["constants"]
    elo = c[2] + c[4]
    return dict(E=d["energy_neutrino"].astype(float), q0=d["q0"].astype(float), Q=c[4],
                w=np.vstack([d[f"w{i}"] for i in range(3)]), conv=c[3] / 1000 * 2 * np.pi, nt=c[5], elo=elo)


def sig_Ebins(s, edges):
    """sigma(E) per channel in E bins, and its MC error (1e-38 cm^2 per nucleus)."""
    idx = np.digitize(s["E"], edges) - 1
    nb = len(edges) - 1
    m = (idx >= 0) & (idx < nb)
    out, err = np.zeros((3, nb)), np.zeros((3, nb))
    thr = s["nt"] * np.diff(edges) / (2000 - s["elo"])
    for k in range(3):
        out[k] = np.bincount(idx[m], weights=s["w"][k][m], minlength=nb) * s["conv"] / thr
        err[k] = np.sqrt(np.bincount(idx[m], weights=s["w"][k][m] ** 2, minlength=nb)) * s["conv"] / thr
    return out, err


def flux_avg(s, lo=150.0):
    m = s["E"] >= lo
    norm = s["conv"] / (s["nt"] * (2000 - lo) / (2000 - s["elo"]))
    v = s["w"][:, m].sum(axis=1) * norm
    e = np.sqrt((s["w"][:, m] ** 2).sum(axis=1)) * norm
    return v, e


summary = {}
edges = np.arange(150, 2001, 100.0)
cen = 0.5 * (edges[1:] + edges[:-1]) / 1000
fig, axes = plt.subplots(3, 4, figsize=(17, 12.5), gridspec_kw=dict(height_ratios=[1, 1, 1.05]))
avail = []
for j, beam in enumerate(["nu", "nubar"]):
    for col_off, nuc in enumerate(["Ar40", "Pb208"]):
        ax = axes[0, 2 * j + col_off]
        a, b = load(f"flat_{nuc}_{beam}_RD0"), load(f"sym_{nuc}_{beam}_RD0")
        if a is None or b is None:
            ax.text(0.5, 0.5, "data missing", transform=ax.transAxes, ha="center"); continue
        sa, ea = sig_Ebins(a, edges); sb, eb = sig_Ebins(b, edges)
        for k in range(3):
            r = sa[k] / sb[k]; re_ = r * np.sqrt((ea[k] / sa[k]) ** 2 + (eb[k] / sb[k]) ** 2)
            ax.plot(cen, r, color=SLOTCOL[k], label=LAB[beam][k])
            ax.fill_between(cen, r - re_, r + re_, color=SLOTCOL[k], alpha=0.2, lw=0)
        rt = sa.sum(0) / sb.sum(0)
        ax.plot(cen, rt, color="k", lw=2, label="total")
        ax.axhline(1, color="grey", lw=0.8)
        ax.set_title(f"{nuc} {'ν' if beam == 'nu' else 'ν̄'}: asymmetric / symmetrised", fontsize=11)
        ax.set_xlabel("E$_\\nu$ [GeV]"); ax.set_ylabel("σ ratio (per channel)")
        ax.legend(fontsize=8, ncol=2); ax.grid(alpha=0.3)
# row 2: flux-averaged ratio per channel and total, all nuclei, nu and nubar
for j, beam in enumerate(["nu", "nubar"]):
    ax = axes[1, 2 * j]; ax2 = axes[1, 2 * j + 1]
    xs = np.arange(len(NUCS))
    for i, nuc in enumerate(NUCS):
        a, b = load(f"flat_{nuc}_{beam}_RD0"), load(f"sym_{nuc}_{beam}_RD0")
        if a is None or b is None:
            continue
        avail.append(f"{nuc}_{beam}")
        va, ea = flux_avg(a); vb, eb = flux_avg(b)
        rs = va / vb; rse = rs * np.sqrt((ea / va) ** 2 + (eb / vb) ** 2)
        rt = va.sum() / vb.sum(); rte = rt * np.sqrt((np.sqrt((ea ** 2).sum()) / va.sum()) ** 2 + (np.sqrt((eb ** 2).sum()) / vb.sum()) ** 2)
        fa, fb = va / va.sum(), vb / vb.sum()
        summary[f"{nuc}_{beam}"] = dict(ratio_slots=rs.tolist(), ratio_slots_err=rse.tolist(), ratio_total=rt, ratio_total_err=rte,
                                        frac_asym=fa.tolist(), frac_sym=fb.tolist(), Q=float(a["Q"]))
        for k in range(3):
            ax.errorbar(i + (k - 1) * 0.18, rs[k], yerr=rse[k], fmt="o", color=SLOTCOL[k], ms=6,
                        label=LAB[beam][k] if i == 0 else None)
        ax.errorbar(i + 0.36, rt, yerr=rte, fmt="s", color="k", ms=7, label="total" if i == 0 else None)
        for k in range(3):
            ax2.bar(i + (k - 1) * 0.27, 100 * (fa[k] - fb[k]), width=0.25, color=SLOTCOL[k],
                    label=LAB[beam][k] if i == 0 else None)
    ax.axhline(1, color="grey", lw=0.8); ax.set_xticks(xs); ax.set_xticklabels([f"{n}\n{(NZ[n][1]-NZ[n][0])/sum(NZ[n]):.2f}" for n in NUCS], fontsize=9)
    ax.set_xlabel("nucleus and (N−Z)/A", fontsize=9)
    ax.set_ylabel("flux-averaged σ ratio, asym / sym"); ax.set_title(f"{'ν' if beam=='nu' else 'ν̄'}: effect of separate p/n Fermi seas (0.15–2 GeV)", fontsize=11)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax2.axhline(0, color="grey", lw=0.8); ax2.set_xticks(xs); ax2.set_xticklabels(NUCS, fontsize=9)
    ax2.set_ylabel("Δ channel fraction [percentage points]"); ax2.set_title(f"{'ν' if beam=='nu' else 'ν̄'}: change of the channel mix (asym − sym)", fontsize=11)
    ax2.legend(fontsize=8); ax2.grid(alpha=0.3)
# row 3: omega shapes (area normalised), asym vs sym, for Ar40 and Pb208, and their ratio
for j, beam in enumerate(["nu", "nubar"]):
    for col_off, nuc in enumerate(["Ar40", "Pb208"]):
        ax = axes[2, 2 * j + col_off]
        a, b = load(f"flat_{nuc}_{beam}_RD0"), load(f"sym_{nuc}_{beam}_RD0")
        if a is None or b is None:
            ax.text(0.5, 0.5, "data missing", transform=ax.transAxes, ha="center"); continue
        oe = np.arange(0, 1.2001, 0.025)
        oc = 0.5 * (oe[1:] + oe[:-1])
        ha = np.histogram((a["q0"] - a["Q"]) / 1000, oe, weights=a["w"].sum(0))[0]
        hb = np.histogram((b["q0"] - b["Q"]) / 1000, oe, weights=b["w"].sum(0))[0]
        ha2 = np.histogram((a["q0"] - a["Q"]) / 1000, oe, weights=a["w"].sum(0) ** 2)[0]
        hb2 = np.histogram((b["q0"] - b["Q"]) / 1000, oe, weights=b["w"].sum(0) ** 2)[0]
        r = (ha / ha.sum()) / (hb / hb.sum())
        re_ = r * np.sqrt(ha2 / ha ** 2 + hb2 / hb ** 2)
        ok = (ha > 0) & (hb > 0)
        ax.plot(oc[ok], r[ok], color="k"); ax.fill_between(oc[ok], (r - re_)[ok], (r + re_)[ok], color="grey", alpha=0.3, lw=0)
        ax.axhline(1, color="grey", lw=0.8); ax.set_ylim(0.6, 1.4)
        ax.set_xlabel("ω = q$_0$ − Q [GeV]"); ax.set_ylabel("shape ratio asym / sym")
        ax.set_title(f"{nuc} {'ν' if beam=='nu' else 'ν̄'}: dσ/dω shape, asym / sym", fontsize=11); ax.grid(alpha=0.3)
fig.suptitle("K. The asymmetric (per-species k$_F$) treatment vs the symmetrised local Fermi gas (old single-k$_F$ physics), same nucleus and same Q", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(f"{FIG}/K_symmetrised.png", dpi=150)
sp = f"{FIG}/summary.json"
s = json.load(open(sp)) if os.path.exists(sp) else {}
s["K_symmetrised"] = summary
s.setdefault("captions", {})["K_symmetrised"] = ("Branch (separate proton and neutron Fermi seas) divided by a control run of the same nucleus with "
    "rho_p = rho_n = rho/2 at every radius (the old single-kF physics), same Q. Isolates the effect of the asymmetric treatment.")
json.dump(s, open(sp, "w"), indent=1)
print("K done; available:", avail)
