#!/usr/bin/env bash
# Run every computation of the paper and report which scripts pass.
# Usage: ./run_all.sh            (the interpreter can be set with PY=/path/to/python)
#        ./run_all.sh --quick    (shorter Section 7 sweep and survey up to n = 20)
set -u
cd "$(dirname "$0")"
PY=${PY:-python3}
QUICK=${1:-}

if [ "$QUICK" = "--quick" ]; then
  jobs=("interval_noncentral_circles.py" "interval_central_circle.py" "interval_gamma_confinement.py"
        "double_point_constants.py" "perturbed_complexes.py" "survey_n_le_50.py 20"
        "t34_second_decomposition.py --quick")
else
  jobs=("interval_noncentral_circles.py" "interval_central_circle.py" "interval_gamma_confinement.py"
        "double_point_constants.py" "perturbed_complexes.py" "survey_n_le_50.py"
        "t34_second_decomposition.py")
fi

status=0
summary=""
for job in "${jobs[@]}"; do
  echo "================================================================ $job"
  start=$(date +%s)
  $PY $job
  rc=$?
  t=$(( $(date +%s) - start ))
  if [ $rc -eq 0 ]; then res="PASS"; else res="FAIL"; status=1; fi
  summary+=$(printf "%-6s %-42s %5ss" "$res" "$job" "$t")$'\n'
done
echo "================================================================ summary"
printf "%s" "$summary"
exit $status
