"""List which expected inputs exist; record in summary.json['inputs_status']."""
import glob
from common import *

E_TOT = [250, 375, 500, 625, 750, 875, 1000, 1250, 1500, 1750, 2000]
exp = {"flat": [], "tot3p3h": [], "tot2p2h": [], "other": []}
for n in NUCS:
    for b in BEAMS:
        exp["flat"].append(f"flat_{n}_{b}_RD0.npz")
for rd in (1, 2):
    for b in BEAMS:
        exp["flat"].append(f"flat_C12_{b}_RD{rd}.npz")
for b in BEAMS:
    for rd in (0, 1, 2):
        for E in E_TOT:
            exp["tot3p3h"].append(f"tot3p3h_C12_{b}_RD{rd}_E{E}.txt")
    for E in [e for e in E_TOT if e <= 1000]:
        exp["tot2p2h"].append(f"tot2p2h_C12_{b}_RD0_E{E}.txt")
for d in ["new_C12_nu", "old_C12_nu", "old_C12_nu_fermioff", "old_C12_nu_fixed", "new_O16_nu", "old_O16_nu", "old_O16_nu_fixed"]:
    exp["other"].append(f"runs/regression/{d}")
exp["other"] += ["harness/mirror_tests.csv", "harness/mirror_tests_old.csv", "data/aux/densities.csv", "data/aux/clampmap.csv", "harness/aux_tools"]
status = {}
for k, lst in exp.items():
    miss = []
    for f in lst:
        if k == "other":
            p = S / f
            ok = (glob.glob(str(p / "events*.root")) != []) if f.startswith("runs/") else p.exists()
        else:
            ok = (DATA / f).exists()
        if not ok:
            miss.append(f)
    status[k] = {"expected": len(lst), "present": len(lst) - len(miss), "missing": miss}
    print(f"{k:8s}: {len(lst) - len(miss):3d}/{len(lst):3d} present" + (f"   missing: {', '.join(miss[:6])}{' ...' if len(miss) > 6 else ''}" if miss else ""))
update_summary({"inputs_status": status})
