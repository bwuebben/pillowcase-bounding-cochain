#!/usr/bin/env bash
# Run all computations of this directory. Usage: ./run_all.sh [--quick]; set PY to choose the interpreter.
set -u
cd "$(dirname "$0")"
status=0
./run_sections_4_5.sh "${1:-}" || status=1
./run_section_7.sh "${1:-}" || status=1
if [ $status -eq 0 ]; then echo "all scripts passed"; else echo "some script failed"; fi
exit $status
