#!/usr/bin/env bash
# Run the computations for Section 7, Remark 6.17 and Appendix A.1-A.3 and report which scripts pass.
# Usage: ./run_section_7.sh            (the interpreter can be set with PY=/path/to/python)
#        ./run_section_7.sh --quick    (shorter runs of kwz_linear_pairings.py and switches_search.py)
# The Khovanov ranks are read from data/khovanov_ranks.txt; to recompute them with Khoca run
#        $PY closures_rational.py --khoca
set -u
cd "$(dirname "$0")"
PY=${PY:-python3}
export PYTHONDONTWRITEBYTECODE=1
QUICK=${1:-}

if [ "$QUICK" = "--quick" ]; then
  linear="kwz_linear_pairings.py --quick"
  search="switches_search.py --quick"
else
  linear="kwz_linear_pairings.py"
  search="switches_search.py"
fi
jobs=("kwz_encoder_check.py" "kwz_hom_check.py" "bn_structure.py" "corner_terms.py" "closures_rational.py"
      "$linear" "sign_positive_t.py" "kwz_traced_curves.py" "$search")

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
