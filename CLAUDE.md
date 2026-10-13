# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

C++17 Monte Carlo event generator / cross-section calculator for neutrino-nucleus 2p2h (and 3p3h) interactions using the Valencia model. Depends on ROOT and GSL. There are no tests and no linter.

## Build and run

```
mkdir -p build && cd build && cmake .. && make -j
```

This produces two executables in `build/`: `program` (the generator, from `main.cc`) and `myroot` (a bare `TRint` interactive ROOT shell that links the `event` dictionaries so saved event files can be read). `clean-all` is a custom clean target.

Run through `./run.sh` **from the repo root**. The script `cd`s into `build/` and expects the already-built `program` there. It also resolves flux files as `../fluxFiles/...`. Its cmake/make lines are commented out, so rebuild by hand. `run.sh` fills in defaults for missing args and accepts the shorthand `t2k` or `minerva` as the 14th argument to select a flux file in `fluxFiles/`.

`program` takes positional arguments (parsed in `main()` in `main.cc`; `run.sh` documents each one):

```
<mode_2p2h> <mode_el> <ReDelta> <emin> <emax> <nr_d> <Global> <target> <efixed> <en_low> <en_high> <nr> [use_flux flux_file flux_hist]
```

- `mode_2p2h`: 0 saves 2p2h events, 1 saves 3p3h events, 2 gives the total 2p2h cross section, 3 gives the total 3p3h cross section.
- `emin`, `emax` and `nr_d` are used only by the total modes. `en_low`, `en_high` and `efixed` are used only by the save modes.
- `target` codes: `0=C12 1=O16 2=O18 3=Na23 4=Ar40 5=Ca40 6=Ca44 7=As75 8=Cd112 9=Pb208`. `Global=1` (global Fermi gas) exists only for C12 and O16.
- `job_eventGen.sh` and `submit_job_*.sub` are Slurm and HTCondor (lxplus) submission templates. Their paths and argument lists are stale or site-specific, so check them against the current argument order before reuse.

## Architecture

- **`main.cc`** parses argv, sets the global physics state, and dispatches to the `saveEvents_*` and `total_*` functions. Each function has two overloads. The argv-driven one is the one `main()` calls. The older interactive one reads from `std::cin`.
- **Global state.** Much of the physics code communicates through globals defined in `lib/constants.h` and similar headers. Examples are `Global`, `CH_2p2h` (CC or CCanti), `ReDelta`, `NR_PART`, `lepton_mass`, `M` and `mpi`. `main()` sets these before dispatch. `main.cc` names its local variable `useGlobal` because a local called `Global` would shadow the global.
- **Local Fermi gas.** `local_fermi_gas(r, nuc)` (`src/utils.cpp`) returns a `FermiGas` struct (`rho_p`, `rho_n`, `rho`, `kF`, `kF_p`, `kF_n`). It is built once per event and passed to the MC and hadron-tensor functions. `kF` is the isospin-averaged value `(3π²ρ/2)^(1/3)` that all current physics uses. `kF_p` and `kF_n` are the per-species values used by the hadron tensors.
- **Isospin slots and species.** The hadron-tensor functions fill three slots `A[0..2]`. For neutrinos they are the final states pp, np, pn (N1', N2'): the W vertex on nucleon N1 emits a pion that is absorbed on N2. For antineutrinos the p<->n mirror gives nn, pn, np. `lib/hadron_tensor_2p2h.h` holds the tables of hole and particle species per slot (`vertex_bubble`) and per interference term (`interference_term`), and `bubble_limits`/`bubble_box` for the hole-energy ranges. Each bubble uses the Fermi momenta of its own species: `DUR` arguments, hole-energy range and Pauli blocking. The binding shift and the Pauli blocking of the direct and crossed nucleon pole follow the species too. The MC samples on the union box over the slots, and the hadron-tensor functions apply each slot's own range and normalisation. The absorbed pion is dressed per charge by `pion_bubbles()` (`PionBubbles`): charged pion (pi+ for neutrinos, pi- for antineutrinos) in a bubble with hole sea kF_hole and particle sea kF_part, pi0 absorbed on a proton or a neutron (so slot 0 and slot 2 get different pi0 numerators), and Delta-hole functions weighted with the N pi -> Delta Clebsch-Gordan coefficients (the crossed term has the weights of the opposite charge). The bubble is `ulind(kF_hole) - crossed(kF_hole) + crossed(kF_part)`, with the new `ulind_crossed`. With `FERMI` on, the charged-pion bubble is evaluated at `w - (kF_hole^2 - kF_part^2)/2M`, which aligns the p and n Fermi energies (for negative shifted energy direct and crossed swap and Im flips); without it the Im part is clamped to <= 0. `FERMI` (hole binding shift, `hadron_tensor_2p2h.cc`) is false by default and `main()` no longer sets it, so every mode runs with it off; set it by hand to turn the shifts on. The Delta width and in-medium Delta propagator still use the total density. For N=Z nuclei (C12, O16, GFG) `pion_bubbles` takes the old `ulind`/`delta_lind` path and the whole output is bit-identical to the old symmetric code.
- **`src/hadron_tensor_2p2h*.cc`** hold the hadron-tensor calculation. The main file is plain 2p2h, `_2b` is the two-body part, and `_rho` is the density-dependent variant. `src/1p1h1pion.cpp` supplies the 1p1h and 1-pion pieces. `lib/MatrixDir.h` and `src/MatrixDir.cc` handle the direct matrix elements (a large header). `lib/Calka.h` is the integration helper.
- **`src/2p2h_MC.cpp` and `src/3p3h_MC.cpp`** are the Monte Carlo integrators and event writers. They cover fixed energy, a uniform energy range, and flux-weighted generation (`saveEvents_2p2h_flux`, which draws Enu from a histogram in a ROOT flux file).
- **Output format.** Events are written to ROOT files through the `event` and `event_3p3h` classes (`lib/event*.h`). CMake generates their dictionaries from `lib/event*LinkDef.h`. `event` stores `kF_p`/`kF_n` next to the averaged `fermi_momentum` (`ClassDef` version 2; old files lack the two branches). Any new branch or field on these classes goes in the header, and `ClassDef` must stay in sync. Output file names encode the options (current, nucleus, energy and the `_ReDelta_*` suffix).
- **Randomness.** `generatormt` is the RNG. `frandom_init(0)` seeds it with `time(NULL)`, so runs are **not** reproducible. To compare two builds, patch `case 0` in `src/generatormt.cpp` to a fixed seed in a scratch copy of the tree, or start the two runs within the same second.

## Nuclei (`lib/utils.h`, `src/utils.cpp`)

Each nucleus is a `nucleus_str`: a density function (`Rho_HO`, `Rho_F` or `Rho_GFG`), its radii and parameters, `Q_2p2h`, and the proton and neutron counts. Targets are looked up by `nucleus_by_code` and `gfg_by_code` and validated by `select_nucleus()` in `main.cc`.

- **Nucleon counts.** `nucleus_str` has `Z` (protons) and `N` (neutrons). The mass number is `Z + N`, which is how `NR_PART` is set.
- **`Q_2p2h`** is S_2p (neutrinos, pp-dominated) and `Q_2p2h_anti` is S_2n (antineutrinos, nn-dominated), AME2020, MeV. `select_nucleus()` copies the nucleus and sets `Q_2p2h` from the current (`CH_2p2h`); np has no own Q. GFG entries keep the legacy value for both. `Q_2P2H_UNSET` is still a sentinel that `select_nucleus()` rejects. Do not use the `Qval` or `Qs` columns for it. They come from the legacy 1p1h table.
- To add a nucleus, add a `nucleus_str` in `src/utils.cpp`, an `extern` in `lib/utils.h`, a case in `nucleus_by_code`, and update `TARGET_CODES` in `main.cc` and the usage text in `run.sh`.

## Conventions

- Energies are in MeV, via the unit constants in `lib/units.h` (e.g. `en *= MeV` after parsing). Densities and radii are in fm units.
- The `nr` argument scales the event count. Its base multiplier differs between modes, so check the specific MC function.
