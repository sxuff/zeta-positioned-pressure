#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-8}"
TABLE_DIR="${TABLE_DIR:-tables}"

# Tables are committed, but regenerating them is the stronger check.
# Set REBUILD_TABLES=1 to force a rebuild even if $TABLE_DIR already exists.
if [[ "${REBUILD_TABLES:-0}" == "1" || ! -f "$TABLE_DIR/w_lower.bin" ]]; then
  python3 src/build_tables.py --output "$TABLE_DIR" --workers "$THREADS"
fi

g++ -O3 -std=c++20 -pthread -ffp-contract=off \
  src/verify_positioned.cpp -o verify_positioned

# Headline certificate.
./verify_positioned 51063 10000000 "$THREADS" "$TABLE_DIR" tangent

# Hardening run: pure interval subdivision, convexity/tangent pruning disabled.
if [[ "${HARDENING:-0}" == "1" ]]; then
  ./verify_positioned 507 100000 "$THREADS" "$TABLE_DIR" no-tangent
fi

python3 src/check_final_bound.py
