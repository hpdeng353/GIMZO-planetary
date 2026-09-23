#!/usr/bin/env bash
# h100_phase6.sh -- end-to-end Phase 6 verification + 0.91+0.11 relaxation on the H100 cluster.
#
# Stages:
#   eostest    compile aneos/test_eos_table.c, run it on the .spheos table, and
#              diff the sample-point output against the Python reference
#              (scripts/test_eos_table.py). Any mismatch = reader bug.
#   build      full planetg build via build_h100.sh (noon1.conf, SYSTYPE=h100)
#   makeplanet submit a SLURM CPU job building the 0.91+0.11 pair, primary 5e6
#   relax      submit one MOONRELAX relaxation SLURM job per built planet
#   all        eostest + build, then makeplanet (relax is submitted by the
#              makeplanet job automatically via --dependency chains)
#
# Usage:  bash scripts/h100_phase6.sh [stage]   (default: all)
# Env overrides: NTASKS (relax MPI ranks, default 64), RELAX_TMAX (default 3.0)

set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$(pwd)
TABLE=$ROOT/impact-out/eos/rock_planet_aneos_62_63.spheos
WORK=$ROOT/runs/091_011
LOG=$WORK/logs
NTASKS=${NTASKS:-48}   # node1: 48 physical cores (96 with HT); nomultithread caps at 48
RELAX_TMAX=${RELAX_TMAX:-3.0}

CONDA_INIT=/home/apps/anaconda3/2024.02/etc/profile.d/conda.sh
py_env() { source "$CONDA_INIT" && conda activate sphexa-planet; }

stage_eostest() {
    echo "=== [eostest] compile C smoke test"
    cc -O2 -o "$WORK/test_eos_table_c" aneos/eos_table.c aneos/test_eos_table.c -lm
    echo "=== [eostest] run C test"
    "$WORK/test_eos_table_c" "$TABLE" | tee "$WORK/eostest_c.out"
    echo "=== [eostest] run Python reference"
    py_env
    python3 scripts/test_eos_table.py "$TABLE" --skip-fnv | tee "$WORK/eostest_py.out"
    echo "=== [eostest] diff sample evaluations (C vs Python)"
    grep -E "^  rho=" "$WORK/eostest_c.out" > "$WORK/eostest_c.samples" || true
    grep -E "^  rho=" "$WORK/eostest_py.out" > "$WORK/eostest_py.samples" || true
    if diff -u "$WORK/eostest_c.samples" "$WORK/eostest_py.samples"; then
        echo "=== [eostest] C and Python sample evaluations AGREE"
    else
        echo "=== [eostest] MISMATCH between C and Python readers -- investigate before continuing" >&2
        exit 1
    fi
    grep -q "RESULT: all checks passed" "$WORK/eostest_c.out"
    grep -q "RESULT: all checks passed" "$WORK/eostest_py.out"
    echo "=== [eostest] PASS"
}

stage_build() {
    echo "=== [build] full planetg build (noon1.conf, SYSTYPE=h100, all cores)"
    py_env   # gsl-config lives in the conda env (GSL_ROOT resolves from it)
    bash build_h100.sh 2>&1 | tee "$WORK/build.log"
    echo "=== [build] PASS: $(ls -la GIZMO)"
}

stage_makeplanet() {
    py_env
    cat > "$WORK/makeplanet.slurm" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=pg-makeplanet
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --hint=nomultithread
#SBATCH --time=08:00:00
#SBATCH --output=$LOG/makeplanet-%j.out
#SBATCH --error=$LOG/makeplanet-%j.err
set -euo pipefail
cd $ROOT
source $CONDA_INIT
conda activate sphexa-planet
python3 scripts/makeplanet_planetg.py \
    --primary-mass-earth 0.91 \
    --secondary-mass-earth 0.11 \
    --primary-particle-count 5000000 \
    --primary-output $WORK/primary_091.h5 \
    --secondary-output $WORK/secondary_011.h5 \
    --force
ls -lh $WORK/primary_091.h5 $WORK/secondary_011.h5
EOF
    local job
    job=$(sbatch --parsable "$WORK/makeplanet.slurm")
    echo "=== [makeplanet] job $job: 0.91 (5e6) + 0.11 M_earth pair"
    echo "$job" > "$WORK/makeplanet.jobid"
}

# relax_one <planet.h5> <name> [dep-jobid]
relax_one() {
    local ic=$1 name=$2 dep=${3:-}
    local rundir=$WORK/relax_$name
    mkdir -p "$rundir"
    # params: copy the noon1 template, drop keys we override, append overrides
    grep -vE "^(InitCondFile|OutputDir|TimeMax|TimeBetSnapshot|EosTable|RelaxTimescale|RelaxUntil|SphericalRelaxUntil|SphericalRelaxReleaseDuration)[[:space:]]" \
        impact-out/init/noon1.params > "$rundir/relax.params"
    cat >> "$rundir/relax.params" <<EOF
InitCondFile                       $ic
OutputDir                          $rundir/
TimeMax                            $RELAX_TMAX
TimeBetSnapshot                    0.1
EosTable                           $TABLE
RelaxTimescale                     0.1
RelaxUntil                         2.0
SphericalRelaxUntil                1.0
SphericalRelaxReleaseDuration      0.5
EOF
    cat > "$rundir/relax.slurm" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=pg-relax-$name
#SBATCH --nodes=1
#SBATCH --ntasks=$NTASKS
#SBATCH --hint=nomultithread
#SBATCH --time=24:00:00
#SBATCH --output=$LOG/relax-$name-%j.out
#SBATCH --error=$LOG/relax-$name-%j.err
set -euo pipefail
module load gcc
module load openmpi
module load hdf5
export OMP_NUM_THREADS=1
# GSL comes from the conda env (no system module on this cluster)
source $CONDA_INIT
conda activate sphexa-planet
export LD_LIBRARY_PATH=\$CONDA_PREFIX/lib:\${LD_LIBRARY_PATH:-}
cd $rundir
mpirun -np $NTASKS $ROOT/GIZMO $rundir/relax.params
EOF
    local depflag=""
    [[ -n "$dep" ]] && depflag="--dependency=afterok:$dep"
    local job
    job=$(sbatch --parsable $depflag "$rundir/relax.slurm")
    echo "=== [relax] $name: job $job (IC=$ic, ntasks=$NTASKS, Tmax=$RELAX_TMAX${dep:+, after $dep})"
}

stage_relax() {
    local dep=""
    [[ -f "$WORK/makeplanet.jobid" ]] && dep=$(cat "$WORK/makeplanet.jobid")
    relax_one "$WORK/primary_091.h5" primary "$dep"
    relax_one "$WORK/secondary_011.h5" secondary "$dep"
    squeue -u "$USER" -o "%.8i %.20j %.8T %.10M %R"
}

mkdir -p "$WORK" "$LOG"
case "${1:-all}" in
    eostest)    stage_eostest ;;
    build)      stage_build ;;
    makeplanet) stage_makeplanet ;;
    relax)      stage_relax ;;
    all)        stage_eostest && stage_build && stage_makeplanet && stage_relax ;;
    *) echo "unknown stage: $1 (eostest|build|makeplanet|relax|all)" >&2; exit 2 ;;
esac
