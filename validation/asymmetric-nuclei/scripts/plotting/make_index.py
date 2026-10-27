"""Write figs/INDEX.md (every figure, caption, key numbers) from summary.json."""
import json
from common import *

sm = json.loads((FIGS / "summary.json").read_text()) if (FIGS / "summary.json").exists() else {}
cap = sm.get("captions", {})
ORDER = ["A_regression", "B_mirror", "C_fig8_nu", "C_fig8_nubar", "D_fractions", "E_diff_C12", "F_neg_C12", "F_neg_pieces",
         "F_neg_nuclei", "G_sigmaA", "G_scaling", "H_isotopes", "H_shapes", "I_radial", "J_clamp"]


def ch(key):
    """final-state channel names (index 0/1/2) for the beam given by the key suffix"""
    return " / ".join(SLOT_FINAL["nubar" if key.endswith("nubar") else "nu"])


def f(x, nd=3):
    return "n/a" if x is None else f"{x:.{nd}g}"


def key_numbers(name):
    L = []
    if name == "A_regression":
        for nuc, d in sm.get("A_regression", {}).items():
            for var, v in d.items():
                r = v["sigma_ratio_new_over_old_pp_np_pn"]
                L.append(f"- {nuc} {var}: sigma ratio new/old pp/np/pn = {f(r[0],5)} / {f(r[1],5)} / {f(r[2],5)}, total {f(v['sigma_ratio_total'],5)}; "
                         f"bit-identical events pp/np/pn = {100*v['identical_fraction_pp_np_pn'][0]:.1f}/{100*v['identical_fraction_pp_np_pn'][1]:.1f}/{100*v['identical_fraction_pp_np_pn'][2]:.1f}%; E_nu aligned: {v['E_aligned']}")
    elif name == "B_mirror":
        b = sm.get("B_mirror", {})
        for k in ["T1_sym_new", "T1_sym_old", "T2_asym", "NC_asym_unmirrored"]:
            if k in b:
                L.append(f"- {k}: max relative difference = {f(b[k]['max_abs_dev'])}, median = {f(b[k]['median_abs_dev'])}")
        if "T1_sym_old" in b:
            L.append(f"- T1_sym_old max relative difference per channel [pp (nu) / nn (nubar), np (nu) / pn (nubar), pn (nu) / np (nubar)] = {[f(x) for x in b['T1_sym_old']['max_dev_per_slot']]} (the pn (nu) / np (nubar) channel lies on the diagonal)")
        for k, v in b.get("T3", {}).items():
            L.append(f"- T3 {k}: median relative difference at eps={[f(e) for e in v['eps']]} = {[f(x) for x in v['median']]}; 90th pct = {[f(x) for x in v['p90']]}; frac(dev>1e-6) = {[f(x) for x in v['frac_dev_gt_1e-6']]}")
    elif name.startswith("C_fig8"):
        beam = "nubar" if name.endswith("nubar") else "nu"
        t = sm.get("C_fig8_table", {}).get(beam, {})
        L.append(f"- paper text band at 1 GeV: {t.get('paper_text_band')}, PRC 102 = {t.get('paper_text_prc102')}")
        for rd in ("RD0", "RD1", "RD2"):
            if rd in t:
                r = t[rd]
                L.append(f"- {rd}: sigma_2p2h(1 GeV) = {f(r.get('sigma_2p2h'))} +- {f(r.get('sigma_2p2h_mc_err'))}, sigma_3p3h = {f(r.get('sigma_3p3h'))}, total = {f(r.get('sigma_total'))}")
        if "our_total_range_RD0_1_2" in t:
            L.append(f"- our total range over RD0/1/2 = {[f(x) for x in t['our_total_range_RD0_1_2']]} vs paper band {t['paper_text_band']}; RD0/paper full curve at 1 GeV = {f(t.get('RD0_over_paper_full_curve_at_1GeV'))}")
        if "flat_over_totalmode_2p2h" in t:
            r = t["flat_over_totalmode_2p2h"]
            L.append(f"- flat-sample 2p2h / total-mode 2p2h at E={r['E']}: {[f(x) for x in r['ratio']]}")
    elif name == "D_fractions":
        for k, v in sm.get("D_fractions", {}).get("flux_avg_fraction", {}).items():
            L.append(f"- {k}: flux-averaged channel fractions ({ch(k)}) = {[round(100*x,1) for x in v]} %")
    elif name == "F_neg_C12":
        for k in ("C12_nu", "C12_nubar"):
            v = sm.get("F_neg_all", {}).get(k)
            if v:
                L.append(f"- {k}: neg. event fraction per channel ({ch(k)}) = {[f(x) for x in v['neg_event_frac_per_slot']]}; negative xsec fraction per channel = {[f(x) for x in v['neg_xsec_frac_per_slot']]}")
    elif name == "F_neg_pieces":
        for k, v in sm.get("F_neg_all", {}).items():
            p = v["piece_neg_over_sigma"]
            L.append(f"- {k}: sum min(w_piece,0)/sum(w): rho {f(p['rho'])}, pi {f(p['pi'])}, int {f(p['int'])}")
    elif name == "F_neg_nuclei":
        for k, v in sm.get("F_neg_all", {}).items():
            L.append(f"- {k}: neg. (event, channel) fraction {f(v['neg_event_frac_pooled'])}, neg. xsec fraction {f(v['neg_xsec_frac_pooled'])}, neg. fraction of events by total weight {f(v['neg_event_frac_total_weight'])}")
    elif name == "G_scaling":
        for k, v in sm.get("G_scaling", {}).get("sigma_per_A_flux_avg", {}).items():
            L.append(f"- {k}: flux-avg sigma/A = {f(v['total'],4)} +- {f(v['err'],2)} (naive pair-luminosity reference {f(v['ref_total'],4)}; ratio {f(v['total']/v['ref_total'],3)})")
        for k, v in sm.get("luminosities", {}).items():
            if isinstance(v, dict):
                L.append(f"- L values {k}: " + ", ".join(f"{a}={f(b,4)}" for a, b in v.items()))
    elif name == "H_isotopes":
        for k, v in sm.get("H_isotopes", {}).get("ratios_flux_avg", {}).items():
            L.append(f"- {k}: flux-avg sigma ratio per channel ({ch(k)}) = {[f(x) for x in v['slots']]}, total {f(v['total'])}; pair-luminosity ratio per channel = {[f(x) for x in v['L_ratio_slots']]}, L_tot ratio {f(v['L_ratio_tot'])}")
    elif name == "H_shapes":
        q = sm.get("H_isotopes", {}).get("Q", {})
        if q:
            L.append("- Q values [MeV]: " + ", ".join(f"{k}={v:.3f}" for k, v in q.items()))
    elif name == "I_radial":
        for k, v in sm.get("I_radial", {}).get("nuclei", {}).items():
            L.append(f"- {k}: kF_n/kF_p centre {f(v.get('kFn_over_kFp_center'))}, max {f(v.get('kFn_over_kFp_max'))}; <r> nu {f(v.get('mean_r_nu'))} fm, nubar {f(v.get('mean_r_nubar'))} fm")
    elif name == "J_clamp":
        for k, v in sm.get("J_clamp", {}).items():
            if isinstance(v, dict) and "frac_events" not in v and "grid_area_fraction_clamped" in v:
                L.append(f"- map {k} (nubar): clamped fraction of the (E_pi,k) grid = {100*v['grid_area_fraction_clamped']:.2f}%")
            if isinstance(v, dict) and "frac_events" in v:
                L.append(f"- {k}: events in clamped region {100*v['frac_events']:.2f}%, |w| in {'pp + np' if k.endswith('_nu') else 'nn + pn'} {100*v['frac_absw_slots01']:.2f}%, |w| in all channels {100*v['frac_absw_all']:.2f}%")
    return L


lines = ["# Validation figures: asymmetric-nuclei 2p2h branch", "",
         "Generated by `scripts/plotting/run_all.sh` (see `SPEC.md`). All cross sections in 1e-38 cm^2 per nucleus unless stated; "
         "flux = flat in E_nu on [E_lo, 2000] MeV; MC errors from sqrt(sum w^2). Machine-readable numbers: `summary.json`.", ""]
n = sm.get("normalisation_check", {})
if n:
    lines += ["## Normalisation check", "",
              f"`{n.get('formula')}` vs the `TOTAL X-sec:` line of err.log: max |relative difference| = {f(n.get('max_abs_rel_diff'))} over {n.get('n_compared')} samples "
              "(differences are at the level of the 6 significant digits printed in err.log).", ""]
    for k, v in n.get("samples", {}).items():
        lines.append(f"- {k}: formula {f(v['sigma_formula'],6)}, err.log {f(v['sigma_errlog'],6)}, rel. diff {f(v['rel_diff'],2)}")
    lines.append("")
ins = sm.get("inputs_status", {})
if ins:
    lines += ["## Input status", ""]
    for k, v in ins.items():
        lines.append(f"- {k}: {v['present']}/{v['expected']} present" + (f"; missing: {', '.join(v['missing'][:12])}{' ...' if len(v['missing']) > 12 else ''}" if v["missing"] else ""))
    lines.append("")
lines += ["## Figures", ""]
for name in ORDER:
    p = FIGS / f"{name}.png"
    lines.append(f"### {name}")
    lines.append("")
    if p.exists():
        lines.append(f"![{name}]({name}.png)")
        lines.append("")
        lines.append(f"**Caption.** {cap.get(name, '(no caption)')}")
        kn = key_numbers(name)
        if kn:
            lines += ["", "Key numbers:"] + kn
    else:
        lines.append("*not produced (inputs missing)*")
    lines.append("")
(FIGS / "INDEX.md").write_text("\n".join(lines))
print("wrote", FIGS / "INDEX.md")
