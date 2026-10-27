#!/bin/bash
# Regenerate every validation figure + summary.json + INDEX.md from whatever inputs exist.
# Usage: run_all.sh            (real data in validation/asymmetric-nuclei/data)
# Missing inputs are skipped gracefully (reported by check_inputs.py and as footnotes on the figures).
cd "$(dirname "$0")" || exit 1
export PYTHONDONTWRITEBYTECODE=1
FIGS="${PLOTS_FIGS:-../figs}"
mkdir -p "$FIGS"
rm -f "$FIGS/summary.json" "$FIGS/INDEX.md"
fail=()
for s in check_inputs check_norm fig_A_beforeafter fig_B_mirror fig_C_nubar fig_D_fractions fig_E_kin fig_F_neg \
         fig_G2_nuclei fig_H_isotopes fig_H_abs fig_I_radial fig_J_clamp fig_K_symmetrised fig_L_C12_Ar40_Ca40 make_index; do
  echo "=================== $s"
  log=$(mktemp)
  nice -n 10 python3 "$s.py" > "$log" 2>&1
  rc=$?
  grep -v "^  \[missing input\]" "$log"
  rm -f "$log"
  if [ "$rc" -ne 0 ]; then echo "!!!! $s FAILED (rc=$rc)"; fail+=("$s"); fi
done
echo "=================== done"
if [ ${#fail[@]} -gt 0 ]; then echo "FAILED scripts: ${fail[*]}"; else echo "all scripts ran"; fi
ls -1 "$FIGS"/*.png 2>/dev/null
