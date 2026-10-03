#!/bin/sh
# Run the computations of Section 6. Pure Python 3, standard library only.
#
#     sh run_all.sh          everything (the q = 11 screens take several minutes)
#     sh run_all.sh fast     the reconstruction checks and the exact q = 7 programs
#
# Exits nonzero if any program reports a failure.
set -e
cd "$(dirname "$0")"
PY=${PY:-python3}

echo "=== reconstruction of Smith's q = 5 curves ==="
for m in tangles resolve earring bigons polygons; do
    echo "--- $m.py"; $PY $m.py
done

echo; echo "=== Alexander polynomials (independent check) ==="
$PY skein_alexander.py

echo; echo "=== exact q = 7 computations in the two-arc algebra ==="
$PY q7_kwz.py --encode
$PY q7_closure_probe.py --slope=-3/4 --selection-certificate
$PY q7_quilt_census.py --all-node-census
$PY q7_quilt_census.py --two-switch-census --strict-pair-census
$PY q7_excluded_orbits.py

echo; echo "=== Smith's perturbed C3 equations, the triple point, the q = 7 trace ==="
$PY c3_perturbed.py
$PY c3_q7_compare.py
$PY surgery_check.py

if [ "$1" = "fast" ]; then
    echo; echo "=== fast part passed ==="; exit 0
fi

echo; echo "=== finite tables ==="
$PY b2_result.py
$PY pretzel_solve.py 3 --triangles-only --max-support 1
$PY pretzel_solve.py 5 --triangles-only --max-support 1
$PY pretzel_solve.py 5 --triangles-only --max-support 1 \
    --blue-epsilon .03 --red-epsilon .10 --red-pinch .25
$PY deform_pent.py
$PY pert_check.py

echo; echo "=== all passed ==="
