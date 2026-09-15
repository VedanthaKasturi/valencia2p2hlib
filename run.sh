#!/bin/bash

# -------------------------------
# Valencia2p2hlib Run Script
# -------------------------------
# Usage:
#   ./run.sh <mode_2p2h> <mode_el> <ReDelta> <emin> <emax> <nr_d> <Global> <target> <efixed> <en_low> <en_high> <nr> [use_flux] [flux_file] [flux_hist]
#
# Description of parameters:
#   mode_2p2h: 0 = save 2p2h, 1 = save 3p3h, 2 = total 2p2h, 3 = total 3p3h
#   mode_el: 0 = CC, 1 = CCanti
#   ReDelta: 0 = 0 MeV, 1 = -50*rho/rho0 MeV, 2 = +30*rho/rho0 MeV
#   emin/emax/nr_d: used in total cross-section modes
#   Global: 0 = Local Fermi Gas, 1 = Global Fermi Gas
#   target: 0=C12 1=O16 2=O18 3=Na23 4=Ar40 5=Ca40 6=Ca44 7=As75 8=Cd112 9=Pb208
#           (Global=1 is only available for 0 and 1)
#   efixed: 1 = fixed Enu (use en_low), 0 = uniform between [en_low, en_high]
#   nr: event-statistics scale factor
#   use_flux: optional (1 to use flux-based generation)
#   flux_file: optional (required if use_flux=1)
#   flux_hist: optional (required if use_flux=1)
#
# Example:
#   ./run.sh 0 0 0 0 0 0 1 0 0 450 650 100 1 t2k
#   ./run.sh 2 0 0 200 2000 100 1 1 0 0 0 0
# -------------------------------

set -e

# Fixed directory for flux files
FLUX_DIR="../fluxFiles"

FLUX_FILE_MINERVA="$FLUX_DIR/MINERvA_ME_MuonNeutrino_flux.root"
FLUX_HIST_MINERVA="flux_hist"
FLUX_FILE_T2K="$FLUX_DIR/t2k_flux.root"
FLUX_HIST_T2K="flux_numu"

# Default arguments (with flux source keyword)
DEFAULT_ARGS=(0 0 0 0 2000 100 0 0 0 450 650 100 1 "t2k")

cd build
#cmake ../
#make -j$(nproc)

# Apply defaults for missing arguments
for i in {1..14}; do
  eval "arg=\${$i}"
  if [[ -z "$arg" ]]; then
    eval "set -- \${@:1:$(($i-1))} \${DEFAULT_ARGS[$(($i-1))]} \${@:$(($i))}"
    echo "Using default for argument $i: ${DEFAULT_ARGS[$(($i-1))]}"
  fi
done

# Interpret shorthand flux keywords ("t2k", "minerva") if only 14 arguments given
if [[ $# -eq 14 ]]; then
  case "${14}" in
    t2k)
      set -- "${@:1:12}" 1 "$FLUX_FILE_T2K" "$FLUX_HIST_T2K"
      ;;
    minerva)
      set -- "${@:1:12}" 1 "$FLUX_FILE_MINERVA" "$FLUX_HIST_MINERVA"
      ;;
  esac
fi

# Validate flux file existence if use_flux is enabled
if [[ "$13" == "1" && ! -f "${14}" ]]; then
  echo "Error: Flux file '${14}' not found."
  exit 1
fi

#if [[ $# -lt 12 ]]; then
#  echo "Error: Not enough arguments provided."
#  echo "Usage: ./run.sh <mode_2p2h> <mode_el> <ReDelta> <emin> <emax> <nr_d> <Global> <target> <efixed> <en_low> <en_high> <nr> [use_flux] [flux_file] [flux_hist]"
#  exit 1
#fi


echo "-----------------------------------"
echo "Launching Valencia2p2hlib generator"
echo "-----------------------------------"
echo "Mode: $1  (0:2p2h, 1:3p3h, 2:tot2p2h, 3:tot3p3h)"
echo "Interaction: $2 (0:CC, 1:CCanti)"
echo "Delta ReShifting: $3"
echo "emin: $4  emax: $5"
echo "nr_d: $6"
echo "Global FG: $7"
echo "Target: $8"
echo "Fixed Enu?: $9"
echo "Enu low: ${10}  high: ${11}"
echo "Event scale (x10^5): ${12}"

if [[ $# -ge 16 ]]; then
  echo "Flux generation enabled"
  echo "Flux file: ${15}"
  echo "Flux hist: ${16}"
fi

echo "-----------------------------------"
./program "$@"
