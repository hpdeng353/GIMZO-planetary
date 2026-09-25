#!/usr/bin/env bash
#SBATCH --job-name=pg-remix
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --hint=nomultithread
#SBATCH --time=12:00:00
#SBATCH --output=/home/hpdeng/planetg/runs/remix/logs/%x-%j.out
#SBATCH --error=/home/hpdeng/planetg/runs/remix/logs/%x-%j.err
# Submit with:  sbatch --export=CASE=square_em tests/remix/run_remix.slurm
# CASE in {square_em, square_es, khi_em, khi_es, rti_em, rti_es}
set -euo pipefail
module load gcc
module load openmpi
module load hdf5
source /home/apps/anaconda3/2024.02/etc/profile.d/conda.sh
conda activate sphexa-planet
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}

CASE=${CASE:?set CASE to square_em|square_es|khi_em|khi_es}
ROOT=/home/hpdeng/planetg
ICS=$ROOT/runs/remix/ics

case "$CASE" in
  square_em) BIN=GIZMO_square; IC=$ICS/square_equal_mass_gizmo;     BASE=square.params ;;
  square_es) BIN=GIZMO_square; IC=$ICS/square_equal_spacing_gizmo;  BASE=square.params ;;
  khi_em)    BIN=GIZMO_khi;    IC=$ICS/khi_iron_rock_gizmo;         BASE=khi.params ;;
  khi_es)    BIN=GIZMO_khi;    IC=$ICS/khi_iron_rock_es_gizmo;      BASE=khi.params ;;
  rti_em)    BIN=GIZMO_rti;    IC=$ICS/rti_iron_rock_em_gizmo;      BASE=rti.params ;;
  rti_es)    BIN=GIZMO_rti;    IC=$ICS/rti_iron_rock_es_gizmo;      BASE=rti.params ;;
  *) echo "unknown CASE $CASE" >&2; exit 2 ;;
esac

d=$ROOT/runs/remix/$CASE
mkdir -p "$d"
sed "s|__ICPATH__|$IC|" "$ROOT/tests/remix/$BASE" > "$d/params.txt"
cd "$d"
mpirun -np 48 "$ROOT/$BIN" params.txt
