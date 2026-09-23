#!/usr/bin/env bash
#SBATCH --job-name=planetg-impact
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --hint=nomultithread
#SBATCH --time=24:00:00
#SBATCH --output=planetg-impact-%j.out
#SBATCH --error=planetg-impact-%j.err

# Submit from the planetg repository root with:
#   sbatch run_impact_h100.slurm
#
# Optional environment overrides, passed with sbatch --export=ALL,NAME=value:
#   PARAMS=init/noon1.params   Parameter file, path relative to RUN_DIR
#   RUN_DIR=impact-out         Self-contained run directory (default: impact-out)
#
# The script records everything needed to reproduce the run inside
# $RUN_DIR/run-info/<date>_<jobid>/: the parameter file, Config.sh, the git
# commit, and the module list. All simulation outputs land in $RUN_DIR itself.

set -euo pipefail

# Initialize environment-modules for non-interactive shells.
if ! type module >/dev/null 2>&1; then
    for module_init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
        if [[ -r "${module_init}" ]]; then
            # shellcheck disable=SC1090
            source "${module_init}"
            break
        fi
    done
fi

module load gcc
module load openmpi
module load hdf5
module load gsl 2>/dev/null || true

repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
run_dir=${RUN_DIR:-impact-out}
params=${PARAMS:-init/noon1.params}
gizmo="${repo_root}/GIZMO"

if [[ ! -x "${gizmo}" ]]; then
    echo "error: ${gizmo} not built; run 'bash build_h100.sh' first" >&2
    exit 1
fi
if [[ ! -f "${repo_root}/${run_dir}/${params}" ]]; then
    echo "error: parameter file '${repo_root}/${run_dir}/${params}' not found" >&2
    exit 1
fi

# --- Record all parameters/code state with the outputs -----------------------
stamp=$(date +%Y%m%d-%H%M%S)
job_tag=${SLURM_JOB_ID:-local}
info_dir="${repo_root}/${run_dir}/run-info/${stamp}_${job_tag}"
mkdir -p "${info_dir}"
cp -f "${repo_root}/${run_dir}/${params}" "${info_dir}/"
cp -f "${repo_root}/Config.sh" "${info_dir}/" 2>/dev/null || true
( cd "${repo_root}" && git rev-parse HEAD ) > "${info_dir}/git-commit.txt" 2>/dev/null || true
module list 2> "${info_dir}/modules.txt" || true
echo "run info recorded in ${info_dir}"

# --- Run ---------------------------------------------------------------------
cd "${repo_root}/${run_dir}"
echo "running ${gizmo} ${params} on ${SLURM_NTASKS:-1} tasks in $(pwd)"
mpirun -np "${SLURM_NTASKS:-1}" "${gizmo}" "${params}"
