# Validation of the `asymmetric-nuclei` branch (26 Oct 2026)

Report: https://claude.ai/artifact/KuBJCsppmZYtGYPP7QT63M (private; `report/report.html` is a local copy of the page).

## plots/ (same order and labels as the report)

| file | made by (scripts/plotting/) |
|---|---|
| 01_A_total_old_vs_joanna | fig_A_beforeafter.py |
| 02_B_charge_symmetry_test | fig_B_mirror.py (data from harness/test_mirror.cc, test_old.cc) |
| 03_C1 … 05_C3 | fig_C_nubar.py |
| 06_D_channel_fractions | fig_D_fractions.py |
| 07_E, 08_E | fig_E_kin.py |
| 09_F, 10_F, 11_F | fig_F_neg.py |
| 12_G1, 13_G1, 14_G2, 15_G3 | fig_G2_nuclei.py |
| 16_H_isotope_ratios, 17_H_isotope_shapes_normalised | fig_H_isotopes.py |
| 18_H_isotope_spectra_absolute | fig_H_abs.py |
| 19_K_asymmetric_vs_symmetrised | fig_K_symmetrised.py |
| 20_I_radial | fig_I_radial.py |
| 21_J_clamp | fig_J_clamp.py (uses harness/aux_tools.cc `clampev`) |
| 22_J_diagram_feynman_bubble, 23_J_diagram_fermi_levels | hand-drawn SVG in the report, rendered to PNG |
| 24_L1_C12_Ar40_Ca40_ownQ, 25_L2_C12_Ar40_Ca40_commonQ | fig_L_C12_Ar40_Ca40.py (numbers in plots/figs/L_C12_Ar40_Ca40.json) |

`run_all.sh` regenerates everything; `common.py` holds shared constants and loaders; `SPEC.md` is the original plotting specification (data formats and normalisation).

## scripts/

- **plotting/** — the figure scripts above.
- **harness/** — standalone C++ programs linked against the branch sources:
  `test_mirror.cc` (charge-symmetry and continuity tests of the hadron tensor, figure B),
  `test_old.cc` (the same test against the old code), `aux_tools.cc` (density profiles, the Im-clamp map and per-event clamp check),
  `dbg_open.cc` (diagnostic for the interference-term steps in B(d)).
- **generation/** — `runjob.sh` runs one generator job and converts its ROOT output with `slim.py`;
  `jobs*.txt` are the exact job lists (tag, seed, beam, program arguments); `patch_test.py` applies the test-only changes.
- **patches/** — every change made to the code copies used for testing (none of it is on any branch):
  - `joanna_branch_instrumented.patch`: RNG seed from `MC_SEED`, optional Q override `Q2P2H_OVERRIDE`, and per-event storage of the ρ/π/interference weight pieces. Weights bit-identical to Joanna's unmodified binary.
  - `joanna_branch_symmetrised_control.patch`: the above plus `SYMMETRIZE_FG` (ρ_p = ρ_n = ρ/2), used only for figure K.
  - `old_code_regression_switches.patch`: seed, Q override, `OLD_FERMI_OFF` and `OLD_INTFIX` switches on the old code, used only to show the regression is bit-identical.

## Data and generator (all inside this folder; nothing lives in scratch)

- **generator/** — git worktree of Joanna's branch at a0e95a6 with `joanna_branch_instrumented.patch` applied; binary `generator/build/program`
  (Release). Rebuild: `cmake -S generator -B generator/build -DCMAKE_BUILD_TYPE=Release && cmake --build generator/build -j8`.
  Checked bit-identical to the 26 Oct production binary (C12 ν, seed 1001, 100k events: E, q0 and all three weights identical).
- **data/** (~3.5 GB) — slimmed event samples `<prefix>_<NUC>_<beam>_RD<k>.npz` (1M events, flat 0–2 GeV, FERMI off) and total-mode `.txt` files:
  `flat_*` (own Q), `sym_*` (symmetrised-gas control, 500k), `qAr_*` (Q2P2H_OVERRIDE = Ar40 Q: 22.757 MeV ν, 16.468 MeV ν̄), `tot*`/`hi*` (total modes), `aux/` (densities, clamp map).
- **runs/** — per-job `cmd.txt` (seed, args, Q override, run time) and logs in `runs/prod/<tag>/`; old-code ROOT files for figures A and C in `runs/beforeafter`, `runs/oldnubar`, `runs/regression`.
- **harness/** — compiled harness programs and their CSV outputs (figure B, J).

Run a job list (each line: `tag seed beam program-args`; `jobs_commonQ.txt` has the Q override as 4th column):

    cd validation/asymmetric-nuclei
    while read tag seed beam args; do bash scripts/generation/runjob.sh $tag $seed $beam $args & done < scripts/generation/jobs.txt; wait
    while read tag seed beam q args; do Q2P2H_OVERRIDE=$q bash scripts/generation/runjob.sh $tag $seed $beam $args & done < scripts/generation/jobs_commonQ.txt; wait

`runjob.sh` (bash) writes to `data/` and `runs/prod/` (override with `DATA_DIR`, `RUNS_DIR`) and uses `generator/build/program` (override with `BIN`, e.g. on a cluster).
`data/`, `runs/`, `harness/` binaries, `generator/` and `plots/figs/` are git-ignored; see `HANDOFF.md` for rebuilding them on a cluster.
Plot scripts read `data/` by default (`PLOTS_DATA` / `PLOTS_S` override) and write working figures to `plots/figs/`; numbered copies go in `plots/`.
