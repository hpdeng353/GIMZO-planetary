#!/usr/bin/env bash
# Rebuild GIZMO_khi / GIZMO_rti with the LONG_Z-fixed confs, then restore the
# production binary. Logs to /tmp/rebuild2.log, done marker /tmp/rebuild2.done
set -euo pipefail
cd /home/hpdeng/planetg
source /home/apps/anaconda3/2024.02/etc/profile.d/conda.sh
conda activate sphexa-planet
CONFIG=conf.khi bash build_h100.sh
mv GIZMO GIZMO_khi
echo "=== GIZMO_khi rebuilt (LONG_Z=0.140625) ==="
CONFIG=conf.rti bash build_h100.sh
mv GIZMO GIZMO_rti
echo "=== GIZMO_rti rebuilt (LONG_Z=0.0703125) ==="
bash build_h100.sh
echo "=== production GIZMO restored ==="
touch /tmp/rebuild2.done
