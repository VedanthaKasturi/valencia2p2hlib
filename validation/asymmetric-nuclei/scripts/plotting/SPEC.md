# Plot spec: validation of the `asymmetric-nuclei` 2p2h branch

All paths below are relative to `S = validation/asymmetric-nuclei/`.
Write all scripts in `S/plots/scripts/`, all figures (PNG, dpi 150, readable fonts >= 9 pt) in `S/plots/figs/`, and a
machine-readable summary of key numbers in `S/plots/figs/summary.json`. Python 3 with numpy, pandas, matplotlib, uproot
(installed). Run python with `python3`. One script per figure group (A..J below) plus a shared `common.py` loader.
Every figure: titled, axes labelled with units, legend, and a one-line caption string saved alongside in
`summary.json["captions"][<figname>]`. Use log scales where values span decades. Use a colour-blind-safe palette and keep
the same colour for the same channel/nucleus across all figures.

## Data

### Flat-flux event samples (main data)
`S/data/flat_<NUC>_<BEAM>_RD<k>.npz`, NUC in {C12,O16,O18,Ar40,Ca40,Ca44,Pb208}, BEAM in {nu,nubar}; RD0 for all,
plus RD1 (Re Sigma_Delta = -50 rho/rho0 MeV) and RD2 (+30) for C12 only. Each has ~1e6 events. Jobs finish over
time: the list of finished jobs is `S/runs/prod/done.txt`. Build and test on what exists; make every script skip
missing files gracefully and print which ones were missing.

Arrays (one entry per stored event; all energies/momenta in MeV, r in fm):
- `energy_neutrino` (E_nu), `q0` (true energy transfer E_nu - E_mu, INCLUDES the Q-value), `qz` (|q|), `theta_lepton`
  (rad), `r`, `rho` (internal units, ignore), `kF_p`, `kF_n`, `E_pi`, `p_pi`, `theta_pi`, `phi_pi`, `E_nuc`, `E_nuc2`, `phi_nuc`.
- weights `w0, w1, w2` for the three isospin slots. Slot meaning by final state (N1', N2'):
  nu: w0 = pp, w1 = np, w2 = pn. nubar: w0 = nn, w1 = pn, w2 = np. (Slot k of nubar is the mirror of slot k of nu.)
  Initial pair: nu slot0 <- (n,p), nu slots 1,2 <- (n,n); nubar slot0 <- (n,p), nubar slots 1,2 <- (p,p).
- pieces `w{i}_rho`, `w{i}_pi`, `w{i}_int` (rho-transverse, pion-longitudinal, interference parts; they sum to `w{i}`).
- `constants` = [nr_events, last E_nu, lepton_mass, C3, Q, n_throws] where C3 = 1e38/cm^2*GeV in internal units and Q is
  the Q-value (MeV) used in that run.

Normalisation (verify it!): sigma in units of 1e-38 cm^2 per nucleus is
`sigma_tot = sum(w) / n_throws * (C3/1000) * 2*pi` (sum over the 3 slots for the total). Check this against the
`TOTAL X-sec:` line printed at the end of `S/runs/prod/<tag>/err.log` (tag = file stem) and record the check in summary.json.
E_nu is uniform on [E_lo, 2000] MeV with E_lo = lepton_mass + Q (= constants[2] + constants[4]); n_throws counts only throws
with E_nu >= E_lo. Hence sigma(E) in an E_nu bin [a,b] (fully above E_lo) is
`sum_{events in bin} w * (C3/1000)*2*pi / (n_throws * (b-a)/(2000-E_lo))`, and the flux-averaged differential
cross section in any variable x is `sum_{bin} w * (C3/1000)*2*pi / (n_throws * dx)`. MC stat error: replace sum(w) by
sqrt(sum(w^2)). Use E_nu bins of 50 MeV starting at the first bin edge above E_lo. Per-nucleon: divide by A
(C12 12, O16 16, O18 18, Ar40 40, Ca40 40, Ca44 44, Pb208 208). Z,N: C12 6,6; O16 8,8; O18 8,10; Ar40 18,22; Ca40 20,20;
Ca44 20,24; Pb208 82,126. Q used: in constants[4].
Define omega = q0 - Q (the energy transfer the hadron tensor sees), Q2 = |q|^2 - q0^2 (MeV^2, show in GeV^2),
E_mu = E_nu - q0, cos(theta_mu) = cos(theta_lepton).

### Total cross sections (other code path)
`S/data/tot3p3h_C12_<BEAM>_RD<k>_E<E>.txt` and `S/data/tot2p2h_C12_<BEAM>_RD0_E<E>.txt`: each line "E sigma" (MeV,
1e-38 cm^2 per nucleus); several lines = independent repeats at the same E -> use mean and standard error of the mean.
3p3h exists for RD0,1,2 and E in {250,375,500,625,750,875,1000,1250,1500,1750,2000}; 2p2h totals for RD0, E<=1000.

### Unit tests (hadron-tensor level)
`S/harness/mirror_tests.csv` and `S/harness/mirror_tests_old.csv`, columns test,point,param,piece,comp,slot,a_re,a_im,b_re,b_im
(read comp as string). Tests: `T1_sym` (new code; a = nu tensor, b = nubar tensor, same slot, symmetric matter, param = kF MeV),
`T1_sym_old` (same test, OLD code), `T2_asym` (a = nu on (kFp,kFn), b = nubar on (kFn,kFp); param = kF_n/kF_p),
`NC_asym_unmirrored` (negative control: b = nubar on the SAME, unmirrored gas), `T3_nu`/`T3_nubar` (a = symmetric path,
b = asymmetric path with relative asymmetry param = eps). pieces: rho, pi, int; comps 00,11,12,03,33; slots 0..2.
For a value use the real part except comp 12 where the imaginary part is the physical one (W3).

### Regression (old branch vs new branch, same seed, same Q)
ROOT files `S/runs/regression/<dir>/events*.root`, tree `treeout`, branches `e/weight_pp`, `e/weight_np`, `e/weight_pn`,
`e/energy_neutrino` (read with uproot `t["e/weight_pp"].array(library="np")`). Dirs: `new_C12_nu` (new branch),
`old_C12_nu` (old branch as delivered: FERMI=true in event mode), `old_C12_nu_fermioff` (old with FERMI off),
`old_C12_nu_fixed` (old with FERMI off + one-line fix of the interference-term normalisation). Same for O16
(`new_O16_nu`, `old_O16_nu`, `old_O16_nu_fixed`). Events are aligned one-to-one (same seed): compare event by event.

### Auxiliary
`S/data/aux/densities.csv` (nucleus,Z,N,r_fm,rho_p_fm3,rho_n_fm3,kF_p_MeV,kF_n_MeV,Q_nu,Q_nubar,rmax).
`S/data/aux/clampmap.csv` (case,kFp,kFn,beam,Epi,kk,ImRaw,ReRaw,ImSym): Im part (MeV^2) of the charged-pion
nucleon-hole bubble before the code clamps Im>0 to 0 (beam 0 = nu pi+, 1 = nubar pi-), ImSym = symmetric-matter
reference. Event-level tool: `S/harness/aux_tools clampev` reads stdin lines "E_pi p_pi kF_p kF_n beam" (MeV; beam 0/1) and
writes "ImRaw ImSym" per line (fast, ~1e6 lines/min is fine; subsample 2e5 events per sample if slow).

### Paper reference: Sobczyk & Nieves, PRC 111, 025502 (2025), Fig. 8 (12C, 2p2h INCLUDING 3p3h, |q|<1.2 GeV)
Exact numbers quoted in the text at E = 1 GeV (1e-38 cm^2): nu band [1.84, 2.61] (+30 and -50 MeV edges),
old model PRC 102 = 1.38; nubar band [0.54, 0.71], old model = 0.31.
Curves read BY EYE from the published figure (uncertainty ~ +-5% or +-0.05, whichever larger; label them
"read by eye from Fig. 8"): (E GeV, sigma)
- nu, full calculation (central, Re Sigma_Delta = 0): (0.25,0.02) (0.375,0.48) (0.5,1.0) (0.625,1.3) (0.75,1.6) (0.875,1.85) (1.0,2.02) (1.25,2.09) (1.5,2.1) (1.75,2.13) (2.0,2.1)
- nu, PRC 102: (0.375,0.2) (0.5,0.5) (0.75,1.0) (1.0,1.38) (1.25,1.47) (1.5,1.48) (2.0,1.49)
- nubar, full calculation: (0.25,0.0) (0.375,0.09) (0.5,0.2) (0.625,0.33) (0.75,0.45) (0.875,0.55) (1.0,0.63) (1.25,0.78) (1.5,0.9) (1.75,0.99) (2.0,1.13)
- nubar, PRC 102: (0.5,0.1) (0.75,0.22) (1.0,0.31) (1.25,0.42) (1.5,0.52) (1.75,0.6) (2.0,0.68)

## Figures

A. `A_regression.png` (C12 and O16 nu): for each of the three old variants vs new, (i) histogram of per-event
   log10|w_new/w_old| for slot pp (or the fraction of bit-identical events written on the panel), (ii) bar chart of
   total sigma ratio new/old per slot (pp, np, pn). Make clear: as delivered (FERMI on) -> FERMI off -> FERMI off + int fix.
B. `B_mirror.png`: 2x2. (a) T1_sym new: scatter of b vs a normalised by max|a| per (point,piece,comp), coloured by slot,
   y=x line; (b) same for T1_sym_old (should show slots 0 and 1 off the diagonal); (c) T2_asym (diagonal) with the
   negative control NC_asym_unmirrored overlaid in grey; (d) T3: median and 90th percentile of |b-a|/max|a| vs eps
   (log-log) for pieces rho, pi (int is ~1e-16; show it too), nu and nubar, with a dashed line proportional to eps.
C. `C_fig8_nu.png`, `C_fig8_nubar.png`: sigma(E) for 12C. Our 2p2h from the flat samples (RD0 central line, RD1/RD2
   band), our 3p3h from the totals (RD0 line, band from RD1/RD2), our 2p2h+3p3h total (interpolate 3p3h to the
   2p2h bin centres linearly), the 2p2h-total-mode points with error bars (independent code path), the paper curves
   read by eye, and the exact text band at 1 GeV as a vertical error bar. Lower panel: our total / paper central.
   Also `C_fig8_table` entries in summary.json: our sigma_2p2h, sigma_3p3h, total at 1 GeV for RD0/1/2 vs paper.
D. `D_fractions.png`: channel fractions vs E_nu: (a) C12 nu (pp,np,pn) and nubar (nn,pn,np) on the same axes with
   mirror channels in the same colour (solid nu, dashed nubar); (b) fraction of the dominant slot (pp for nu, nn for
   nubar) vs E_nu for all nuclei; (c) bar chart of flux-averaged fractions per nucleus and beam.
E. `E_diff_C12.png`: C12 flux-averaged per-slot differential cross sections, nu (solid) vs nubar (dashed),
   mirror slots same colour: dsigma/dq0, dsigma/domega, dsigma/d|q|, dsigma/dQ2, dsigma/dE_mu, dsigma/dcos(theta_mu).
F. Negative weights: `F_neg_C12.png`: for C12 nu and nubar, per slot: (a) fraction of events with w<0 vs E_nu,
   (b) negative fraction of the cross section, sum_{w<0}|w| / sum(w), vs E_nu, (c) same vs q0, (d) vs |q|.
   `F_neg_pieces.png`: which piece is negative: for each sample (all nuclei, both beams) the sum of negative
   contributions of each piece relative to the total sigma: sum(min(w_piece,0))/sum(w), stacked bars; and the
   distribution of w_int/w for events with w<0. `F_neg_nuclei.png`: negative event fraction and negative xsec
   fraction per nucleus and beam vs (N-Z)/A, and vs r and vs local kF_n/kF_p for Ar40 and Pb208.
G. Nuclei: `G_sigmaA.png`: sigma/A vs E_nu for all nuclei (nu and nubar panels), and ratio to C12;
   `G_scaling.png`: flux-averaged sigma/A vs A for nu and nubar with, as reference, the naive pair-luminosity
   expectation normalised to C12: L/A with L = 4*pi*int r^2 rho_a(r) rho_b(r) dr from densities.csv (total:
   L_tot = int rho^2 with rho = rho_p + rho_n); also per slot: nu slot0 vs L_np, nu slots1+2 vs L_nn, nubar slot0 vs
   L_np, nubar slots1+2 vs L_pp. Write the L values into summary.json. Label the reference "naive rho_a*rho_b pair
   counting (no Pauli blocking), normalised to C12".
H. Asymmetric isotopes: `H_isotopes.png`: per-slot sigma ratios O18/O16, Ca44/Ca40 and Ar40/Ca40 vs E_nu, for nu and
   nubar, with the pair-luminosity ratio for that slot as a horizontal dashed line. `H_shapes.png`: area-normalised
   dsigma/dq0 and dsigma/domega for O16 vs O18 and Ca40 vs Ca44 (nu and nubar), with the Q values written in the legend,
   to separate the rigid Q shift from hadron-tensor changes.
I. Radial: `I_radial.png`: per nucleus, (a) density profiles rho_p(r), rho_n(r) and kF_n/kF_p(r) (twin axis),
   (b) dsigma/dr normalised (nu and nubar), (c) slot fractions vs r for Ar40, Ca40, Pb208 (nu and nubar).
J. Clamp: `J_clamp.png`: (a) map of ImRaw over (E_pi, k_pi) for nubar at the Pb208 and Ar40 centres with the
   ImRaw>0 (clamped) region highlighted; (b) for every flat sample, fraction of events (and of |w| in slots 0,1)
   whose (E_pi,p_pi,kF_p,kF_n) lie in the clamped region (use clampev on a subsample).

Finally write `S/plots/figs/INDEX.md` listing every figure with its caption and the key numbers.
