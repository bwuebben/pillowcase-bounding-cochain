#!/usr/bin/env bash
# Run the computations for Sections 4 and 5 and Appendix A.3 and report which scripts pass.
# Usage: ./run_sections_4_5.sh            (the interpreter can be set with PY=/path/to/python)
#        ./run_sections_4_5.sh --quick    (shorter run of sheared_numerics.py)
# pretzel_chirality_snappy.py needs SnapPy; it is reported as SKIP if SnapPy is not installed.
set -u
cd "$(dirname "$0")"
PY=${PY:-python3}
export PYTHONDONTWRITEBYTECODE=1
QUICK=${1:-}

if [ "$QUICK" = "--quick" ]; then
  numerics="sheared_numerics.py --quick"
else
  numerics="sheared_numerics.py"
fi
jobs=("straightness_identity.py" "sheared_polygon_counts.py" "pretzel_groups.py" "twisted_alexander.py"
      "sheared_shear_check.py" "sheared_resolution_type.py" "$numerics" "sheared_split.py"
      "twisted_numerics.py" "pretzel_chirality_snappy.py")

status=0
summary=""
for job in "${jobs[@]}"; do
  echo "================================================================ $job"
  start=$(date +%s)
  $PY $job
  rc=$?
  t=$(( $(date +%s) - start ))
  if [ $rc -eq 0 ]; then res="PASS"; elif [ $rc -eq 2 ]; then res="SKIP"; else res="FAIL"; status=1; fi
  summary+=$(printf "%-6s %-42s %5ss" "$res" "$job" "$t")$'\n'
done
echo "================================================================ summary"
printf "%s" "$summary"
exit $status
