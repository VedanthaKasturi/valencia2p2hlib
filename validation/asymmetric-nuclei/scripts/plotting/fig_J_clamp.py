"""Figure J: the Im>0 clamp of the charged-pion bubble: maps and the fraction of events/|w| sitting in the clamped region."""
import subprocess
import zlib
import pandas as pd
from matplotlib.colors import SymLogNorm
from common import *

NSUB = 500_000
TOOL = HARNESS / "aux_tools"


def clamp_flags(s):
    """ImRaw, ImSym and subsample index for a Sample (cached in plots/cache)."""
    cf = CACHE / f"clamp_{s.tag}_{os.environ.get('PLOTS_DATA_TAG', 'real')}.npz"
    if cf.exists() and cf.stat().st_mtime > s.path.stat().st_mtime and np.load(cf)["nsub"] == min(NSUB, s.n):
        z = np.load(cf)
        return z["idx"], z["raw"], z["sym"]
    rng = np.random.default_rng(zlib.crc32(s.tag.encode()))
    idx = np.sort(rng.choice(s.n, min(NSUB, s.n), replace=False))
    beam = 0 if s.beam == "nu" else 1
    cols = np.column_stack([s.a("E_pi")[idx], s.a("p_pi")[idx], s.a("kF_p")[idx], s.a("kF_n")[idx], np.full(len(idx), beam)])
    inp = "\n".join(" ".join(f"{v:.9g}" for v in row) for row in cols) + "\n"
    r = subprocess.run(["nice", str(TOOL), "clampev"], input=inp, capture_output=True, text=True, check=True)
    v = np.array([l.split() for l in r.stdout.strip().split("\n")], dtype=float)
    assert v.shape == (len(idx), 2), v.shape
    np.savez(cf, idx=idx, raw=v[:, 0], sym=v[:, 1], nsub=len(idx))
    return idx, v[:, 0], v[:, 1]


# ------------------------------------------------------------------------------------------ (b) events
res = {}
for nuc in NUCS:
    for beam in BEAMS:
        s = load_sample(nuc, beam)
        if s is None:
            continue
        if not TOOL.exists():
            note_missing("harness/aux_tools")
            break
        idx, raw, sym = clamp_flags(s)
        cl = raw > 0
        w = [np.abs(s.w(k))[idx] for k in range(3)]
        w01 = w[0] + w[1]; wall = w01 + w[2]
        # also the signed weight in the clamped region
        res[(nuc, beam)] = {"n_sub": len(idx), "frac_events": float(cl.mean()),
                            "frac_absw_slots01": float(w01[cl].sum() / w01.sum()),
                            "frac_absw_all": float(wall[cl].sum() / wall.sum()),
                            "frac_absw_per_slot": [float(w[k][cl].sum() / w[k].sum()) for k in range(3)],
                            "frac_events_symmetric_ref_Im_gt_0": float((sym > 0).mean()),
                            "kF_n_over_kF_p_median": float(np.median(s.a("kF_n")[idx] / np.maximum(s.a("kF_p")[idx], 1e-9)))}
        s.free("E_pi", "p_pi", "kF_p", "kF_n", "w0", "w1", "w2", "wtot")
        del s
miss = missing_list()
update_summary({"J_clamp": {f"{n}_{b}": v for (n, b), v in res.items()},
                "J_clamp_note": f"clamped = raw Im of the charged-pion bubble > 0 (set to 0 by the code); evaluated with aux_tools clampev on a random subsample of {NSUB} events per sample"})

# ------------------------------------------------------------------------------------------ figure
fig = plt.figure(figsize=(16, 13))
gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1], width_ratios=[1, 1, 0.04], hspace=0.28, wspace=0.18)
cm_path = AUX / "clampmap.csv"
cmap_df = pd.read_csv(cm_path) if cm_path.exists() else None
if cmap_df is None:
    note_missing("aux/clampmap.csv")
extra = {}
last_pcm = None
for j, case in enumerate(["Pb208_r0", "Ar40_r0"]):
    ax = fig.add_subplot(gs[0, j])
    if cmap_df is None:
        no_data(ax); continue
    d = cmap_df[(cmap_df.case == case) & (cmap_df.beam == 1)]
    piv = d.pivot_table(index="kk", columns="Epi", values="ImRaw")
    E_, K_ = np.meshgrid(piv.columns.to_numpy(), piv.index.to_numpy())
    Z = piv.to_numpy()
    norm = SymLogNorm(linthresh=1e3, vmin=-np.nanmax(np.abs(Z)), vmax=np.nanmax(np.abs(Z)), base=10)
    last_pcm = ax.pcolormesh(E_, K_, Z, cmap="RdBu_r", norm=norm, shading="auto", rasterized=True)
    mask = (Z > 0).astype(float)
    ax.contourf(E_, K_, mask, levels=[0.5, 1.5], colors="none", hatches=["////"])
    ax.contour(E_, K_, mask, levels=[0.5], colors="k", linewidths=1.2)
    kfp, kfn = d.kFp.iloc[0], d.kFn.iloc[0]
    ax.set_xscale("log"); ax.set_xlim(0.5, 400)
    ax.set_xlabel(r"$E_\pi$ [MeV]  (log axis)"); ax.set_ylabel(r"$k_\pi$ [MeV]")
    ax.set_title(f"(a) {case.split('_')[0]}, r=0:  $k_{{F,p}}$={kfp:.0f}, $k_{{F,n}}$={kfn:.0f} MeV;  $\\bar\\nu$ ($\\pi^-$)", fontsize=10.5)
    area = float((Z > 0).mean())
    ax.text(0.02, 0.97, f"hatched: Im > 0 (clamped), {100*area:.0f}% of the grid", transform=ax.transAxes, va="top", fontsize=9.5,
            bbox=dict(fc="white", ec="none", alpha=0.85))
    extra[case] = {"grid_area_fraction_clamped": area, "kFp": float(kfp), "kFn": float(kfn)}
    # overlay event density of the matching sample (E_pi, p_pi), if present
    nuc = case.split("_")[0]
cax = fig.add_subplot(gs[0, 2])
if last_pcm is not None:
    cb = fig.colorbar(last_pcm, cax=cax)
    cb.set_label(r"Im (raw) of the $\pi N^{-1}$ bubble [MeV$^2$] (symlog)")
else:
    cax.axis("off")
# (b) -------------------------------------------------------------------------------------------------------------------
axb = fig.add_subplot(gs[1, :2])
if res:
    labels, x = [], 0
    wd = 0.26
    for beam in BEAMS:
        for nuc in NUCS:
            if (nuc, beam) in res:
                r = res[(nuc, beam)]
                axb.bar(x - wd, 100 * r["frac_events"], wd, color="#555555", label="events" if x == 0 else None)
                axb.bar(x, 100 * r["frac_absw_slots01"], wd, color=SLOT_COL[0], label=r"$|w|$ in pp + np ($\nu$) / nn + pn ($\bar\nu$)" if x == 0 else None)
                axb.bar(x + wd, 100 * r["frac_absw_all"], wd, color="#E69F00", label=r"$|w|$ in all channels" if x == 0 else None)
                for dx, key in ((-wd, "frac_events"), (0, "frac_absw_slots01"), (wd, "frac_absw_all")):
                    axb.text(x + dx, 100 * r[key], f"{100 * r[key]:.2g}", ha="center", va="bottom", fontsize=9, rotation=90)
            else:
                axb.text(x, 0.5, "n/a", rotation=90, ha="center", va="bottom", color="#a00000", fontsize=9)
            labels.append(f"{nuc}\n{'nu' if beam == 'nu' else 'nubar'}")
            x += 1
        x += 0.6
        labels.append("")
    xs = []
    xx = 0
    for beam in BEAMS:
        for nuc in NUCS:
            xs.append(xx); xx += 1
        xs.append(np.nan); xx += 0.6
    tk = [(a, l) for a, l in zip(xs, labels) if np.isfinite(a)]
    axb.set_xticks([a for a, _ in tk]); axb.set_xticklabels([l for _, l in tk], fontsize=9)
    axb.set_ylabel("share inside the clamped region  [%]")
    axb.set_title(f"(b) fraction of events and of $|w|$ whose ($E_\\pi$, $p_\\pi$, $k_{{F,p}}$, $k_{{F,n}}$) lie where the raw Im of the charged-pion bubble is > 0  "
                  f"(subsample of {NSUB:,} events)", fontsize=10.5)
    axb.legend(loc="upper left", fontsize=9.5)
    axb.margins(y=0.15)
    axb.grid(axis="x", visible=False)
else:
    no_data(axb, "no flat samples available yet")
stamp_missing(fig, miss)
fig.suptitle(r"J. The Im$>$0 clamp in the charged-pion nucleon-hole bubble (asymmetric gas)", fontsize=12.5, y=0.985)
fig.subplots_adjust(left=0.07, right=0.93, top=0.93, bottom=0.08)
save_fig(fig, "J_clamp", "Clamp of the charged-pion nucleon-hole bubble (Im>0 set to 0): (a) raw Im over (E_pi, k_pi) for nubar at the Pb208 and Ar40 centres with the clamped region hatched; (b) share of events and of |w| (pp + np for nu, nn + pn for nubar; and all channels) of every flat sample inside the clamped region.", extra)
finish("J")
