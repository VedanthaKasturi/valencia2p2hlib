# Convert one 2p2h event ROOT file into a compact npz for analysis. Usage: slim.py file.root out.npz beam(nu|nubar)
import sys, numpy as np, uproot
f = uproot.open(sys.argv[1]); t = f["treeout"]; beam = sys.argv[3]
c = f["constants"].member("fElements")
wk = ["weight_pp","weight_np","weight_pn"] if beam=="nu" else ["weight_anti_nn","weight_anti_pn","weight_anti_np"]
kin = ["energy_neutrino","q0","qz","theta_lepton","r","rho","kF_p","kF_n","E_pi","p_pi","theta_pi","phi_pi","E_nuc","E_nuc2","phi_nuc"]
out = {"constants": np.asarray(c, dtype=np.float64)}
for k in kin: out[k] = t[f"e/{k}"].array(library="np").astype(np.float32)
for i,k in enumerate(wk): out[f"w{i}"] = t[f"e/{k}"].array(library="np").astype(np.float64)
for p in ["rho","pi","int"]:
    a = t[f"e/w_{p}[3]"].array(library="np")
    for i in range(3): out[f"w{i}_{p}"] = a[:,i].astype(np.float64)
np.savez_compressed(sys.argv[2], **out)
print("slimmed", sys.argv[2], len(out["w0"]))
