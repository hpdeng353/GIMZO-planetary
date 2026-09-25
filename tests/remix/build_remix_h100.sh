#!/usr/bin/env bash
# Build the two REMIX test binaries (square: ideal gas; khi: EOS_ANEOS) and
# restore the production build afterwards. ICs are converted locally with
# scripts/remix_ic_to_gizmo.py and uploaded to runs/remix/ics/ beforehand.
#
# Usage: bash tests/remix/build_remix_h100.sh
set -euo pipefail
cd "$(dirname "$0")/../.."

CONFIG=conf.square bash build_h100.sh
mv GIZMO GIZMO_square
echo "=== GIZMO_square built ==="

CONFIG=conf.khi bash build_h100.sh
mv GIZMO GIZMO_khi
echo "=== GIZMO_khi built ==="

CONFIG=conf.rti bash build_h100.sh
mv GIZMO GIZMO_rti
echo "=== GIZMO_rti built ==="

bash build_h100.sh   # restore production Config.sh (moon.conf) + GIZMO
echo "=== production GIZMO restored ==="
