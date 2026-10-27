# Handoff: asymmetric-nuclei 2p2h validation → cluster production

Written 2026-10-27 for the Claude Code session that continues this work on the cluster. Read this first, then `README.md`.

## Where we are

- `feature/NewNuclei` now contains Joanna Sobczyk's `asymmetric-nuclei` branch (fast-forward to a0e95a6, 27 Oct 2026):
  per-species Fermi momenta kF_p/kF_n, charge-dependent pion dressing, the antineutrino fix, S_2p (ν) / S_2n (ν̄) as Q_2p2h,
  kF_p/kF_n stored in the event record. See the root `CLAUDE.md` for the code layout.
- A validation of that branch was done on a laptop (26–27 Oct) with 1M events per nucleus and beam. All plots are in `plots/`
  (numbered 01–25, same labels as the report). **The goal on the cluster is to redo these plots with more statistics.**
- The bulk data (3.5 GB of `.npz` samples, ROOT files for figures A and C, harness CSVs) stayed on the laptop and is git-ignored.
  Regenerate it on the cluster from the job lists below.

## Physics setup used for every sample

- Branch code at a0e95a6 plus the **test-only** patch `scripts/patches/joanna_branch_instrumented.patch`
  (RNG seed from env `MC_SEED`, Q override from env `Q2P2H_OVERRIDE` in MeV, per-event ρ/π/interference weight pieces `w_rho/w_pi/w_int`).
  The patch does not change any weight (checked bit-identical to the unpatched binary). Never commit it to a branch.
- FERMI off (branch default), Re Σ_Δ = 0 (`ReDelta` arg 0) unless the tag says RD1/RD2, flat E_ν 0–2 GeV, event mode (`mode_2p2h 0`).
- "Own Q" = the Q values the branch ships (`src/utils.cpp`). Those are the physics prediction.
  `qAr_*` runs put C12 and Ca40 at the Ar40 Q (22.757 MeV ν, 16.468 MeV ν̄) to separate the Q effect from nuclear structure.
  They are a diagnostic only, not a prediction.

## Rebuilding the samples on the cluster

1. Build the instrumented generator in a separate worktree (keeps the branch clean):

       git worktree add --detach validation/asymmetric-nuclei/generator a0e95a6
       cd validation/asymmetric-nuclei/generator && git apply ../scripts/patches/joanna_branch_instrumented.patch
       cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j

   Needs ROOT and GSL. `slim.py` needs Python 3 with numpy and uproot.
2. Each job-list line is `tag seed beam program-args` (`jobs_commonQ.txt` has the Q override as the 4th column).
   One line = one call of `scripts/generation/runjob.sh` (bash). It runs the generator in `runs/prod/<tag>/`,
   slims `events*.root` into `data/<tag>.npz` and deletes the ROOT file. Override `DATA_DIR`, `RUNS_DIR`, `BIN` as needed.
   - `jobs.txt`: production (`flat_*` for 7 nuclei × ν/ν̄, RD0/1/2 for C12, total-mode `tot*`).
   - `jobs_sym.txt`: symmetrised-gas control for figure K (needs `joanna_branch_symmetrised_control.patch` instead).
   - `jobs_hi.txt`, `jobs_q.txt`: high-precision total-mode points.
   - `jobs_commonQ.txt`: `qAr_*` runs for figure L2.
3. More statistics: the last program argument `nr` scales the event count (10 = 1M events in these modes).
   Better to run several jobs per tag with **different seeds** and concatenate, than one huge job. The plotting
   loaders read a single `<tag>.npz`, so merge the pieces before plotting: concatenate the arrays, keep `constants`,
   and **sum** `constants[0]` (nr_events) and `constants[5]` (n_throws, which the σ normalisation divides by).
   Then run `check_norm.py`.
   The submission templates in the repo root (`submit_job_*.sub`, `job_eventGen.sh`) are stale or site-specific.
4. Timing on the laptop with 1M events: C12 ≈ 6–12 min, Ar40/Ca40 ≈ 9–17 min per job.
5. Reproducibility check before trusting a cluster build:
   `bash scripts/generation/runjob.sh check_C12 1001 nu 0 0 0 0 0 0 0 0 0 0 2000 1` (C12 ν, 100k events). The laptop build gives
   first three `energy_neutrino` = 706.8912, 445.3259, 524.8417 MeV and weight sums
   Σw0 (pp) = 5.920225562498783e-12, Σw1 (np) = 4.184293914062319e-13, Σw2 (pn) = 1.0558053841334325e-12.
   A different compiler or libm can break bit-equality; then the sums should agree to well below 0.1 %.

## Plotting

`scripts/plotting/run_all.sh` regenerates every figure from `data/` into `plots/figs/` (defaults resolve relative to this folder;
`PLOTS_S`, `PLOTS_DATA`, `PLOTS_FIGS` override). Copy the finals into `plots/` with the numbered names (table in `README.md`).
Figures A and C also need old-code ROOT runs (`runs/beforeafter`, `runs/oldnubar`) and figure B/J need harness outputs.
Those are regression checks against the old code and do not need more statistics. Skip them unless asked.

## Reference results to reproduce (laptop, 1M events)

- Flux-averaged (0.15–2 GeV) scaling: σ_ν ∝ A^1.11, N^1.05; σ_ν̄ ∝ A^1.01, Z^1.08.
- σ/A ratio to C12, own Q: ν Ar40 1.130 ± 0.015, Ca40 1.123 ± 0.014; ν̄ Ar40 1.049 ± 0.016, Ca40 1.073 ± 0.016.
  At the common Ar40 Q: ν Ar40/Ca40 = 1.02 ± 0.02; ν̄ Ar40/Ca40 = 0.89 ± 0.02.
- Separate p/n Fermi seas vs a symmetrised gas (same Q): Ar40 +7 % ν / −7 % ν̄; Pb208 +13 % / −21 %.
- Negative weights cancel 2.5–3.3 % of σ for all nuclei. C12 ν pp share 80 %; C12 ν̄ nn 77 % / pn 7.5 % / np 16 %.

## Literature to compare with

- Martinez-Consentino, Segovia, Amaro, arXiv:2509.08916: per-nucleon 2p2h MEC response Ar40/C12 ≈ 1.04–1.30 (RMF), 1.33–1.53 (RFG);
  Ca40 as an Ar40 proxy underestimates ν 2p2h by about 10 %. Ours gives about 2 % for ν. This is open and worth raising with Joanna.
- Amaro et al., PRC 95, 065502 (2017): 2p2h peak ∝ A·k_F². Colle et al., PRC 92, 024604 (2015): SRC pairs ∝ A^1.12.
  Mosel & Gallmeister, PRC 94, 034610 (2016): 2N knockout follows ⟨ρ⟩ = 0.145 − 0.147 A^(−1/3).
- No literature check exists yet for the ν̄ A-scaling.

## Open points for Joanna (from the validation)

Q per channel (np has no own Q); interference terms `tr3_mu tr2_nu` / `tr2_mu tr3_nu` seem added to the swapped np/pn entries
(inherited, < 1 % of the total); no Fermi-level alignment with FERMI off; ask for the numerical PRC 111 Fig. 8 curves;
total-σ modes reuse `nr_d` and overwrite it above 1 GeV.

## How Vedantha wants results presented

- Answer exactly the question, briefly. Physics labels pp/np/pn (ν) and nn/pn/np (ν̄), never internal names like "slot".
- σ(E) with a ratio panel. Minimal comparisons, only the ones asked for. No 3p3h added. No naive pair-counting curves.
- Prefer plots over tables. Store final PNGs numbered in `plots/`.
- Keep every output inside the project (never in a temporary or scratch directory). Record the seed, args and Q override next to each output (`runjob.sh` does this in `cmd.txt`).
- Ask before running anything beyond what was requested (e.g. extra diagnostic runs).
