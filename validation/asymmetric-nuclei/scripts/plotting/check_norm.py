"""Normalisation check: sigma_tot = sum(w)/n_throws*(C3/1000)*2*pi  vs  'TOTAL X-sec:' in err.log (SPEC, Data)."""
import glob
import re
import sys
import numpy as np
from common import *


def logged_total(path):
    if not path.exists():
        return None
    m = re.findall(r"TOTAL X-sec:\s*([0-9.eE+-]+)", path.read_text())
    return float(m[-1]) if m else None


rows = {}
# 1. finished flat jobs (npz + err.log)
for nuc in NUCS:
    for beam in BEAMS:
        for rd in (0, 1, 2):
            tag = f"flat_{nuc}_{beam}_RD{rd}"
            if rd > 0 and nuc != "C12":
                continue
            p = DATA / f"{tag}.npz"
            if not p.exists():
                note_missing(p.name)
                continue
            s = Sample(nuc, beam, rd)
            lg = logged_total(RUNS / "prod" / tag / "err.log")
            mine = s.sigma_tot()
            rows[tag] = {"sigma_formula": mine, "sigma_errlog": lg,
                         "rel_diff": None if not lg else mine / lg - 1.0,
                         "per_slot": [s.sigma_tot(k) for k in range(3)], "n_events": s.n,
                         "n_throws": s.n_throws, "Q": s.Q, "E_lo": s.E_lo}
# 2. the regression 'new' event runs (same program, 1e5 events) -- available from the start
import uproot
for d in sorted(glob.glob(str(RUNS / "regression" / "new_*_nu*"))):
    tag = Path(d).name
    fn = glob.glob(d + "/events*.root")
    if not fn:
        continue
    f = uproot.open(fn[0])
    t = f["treeout"]
    c = np.asarray(f["constants"].member("fElements"), float)
    w = [t[f"e/{k}"].array(library="np") for k in ["weight_pp", "weight_np", "weight_pn"]]
    mine = sum(x.sum() for x in w) / c[5] * (c[3] / 1000) * 2 * np.pi
    lg = logged_total(Path(d) / "err.log")
    rows["regression/" + tag] = {"sigma_formula": float(mine), "sigma_errlog": lg,
                                 "rel_diff": None if not lg else float(mine / lg - 1.0),
                                 "n_events": len(w[0]), "n_throws": c[5], "Q": c[4]}
ok = [r["rel_diff"] for r in rows.values() if r["rel_diff"] is not None]
summ = {"normalisation_check": {
    "formula": "sigma_tot = sum(w0+w1+w2)/n_throws*(C3/1000)*2*pi  [1e-38 cm^2 per nucleus]",
    "samples": rows,
    "max_abs_rel_diff": float(np.max(np.abs(ok))) if ok else None,
    "n_compared": len(ok)}}
update_summary(summ)
print(f"{'sample':38s} {'formula':>10s} {'err.log':>10s} {'rel.diff':>10s}")
for k, r in rows.items():
    print(f"{k:38s} {r['sigma_formula']:10.5f} {r['sigma_errlog'] if r['sigma_errlog'] is not None else float('nan'):10.5f} "
          f"{r['rel_diff'] if r['rel_diff'] is not None else float('nan'):10.2e}")
print("max |rel diff| =", summ["normalisation_check"]["max_abs_rel_diff"], "over", len(ok), "samples")
