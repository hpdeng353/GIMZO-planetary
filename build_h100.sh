#!/usr/bin/env bash
# Build planetg (GIZMO) on the H100 cluster using all available CPU cores.
#
# Usage:
#   bash build_h100.sh                 # build with impact.conf
#   CONFIG=conf.khi bash build_h100.sh  # build another tracked configuration
#   JOBS=32 bash build_h100.sh         # override parallel jobs (default: nproc)
#
# The script loads the GNU/OpenMPI/HDF5(/GSL) modules, selects SYSTYPE="h100"
# in the Makefile, copies the chosen .conf to Config.sh, and rebuilds from
# clean (GIZMO requires a full rebuild after any Config.sh change).

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
module load gsl 2>/dev/null || true   # gsl may not exist as a module

# --- Resolve HDF5_ROOT -------------------------------------------------------
: "${HDF5_ROOT:=${HDF5_DIR:-${HDF5_HOME:-}}}"
if [[ -z "${HDF5_ROOT}" ]]; then
    for helper in h5pcc h5cc; do
        if command -v "${helper}" >/dev/null 2>&1; then
            HDF5_ROOT=$(cd -- "$(dirname -- "$(command -v "${helper}")")/.." && pwd)
            break
        fi
    done
fi
if [[ -z "${HDF5_ROOT}" || ! -d "${HDF5_ROOT}/include" ]]; then
    echo "error: cannot locate HDF5; set HDF5_ROOT explicitly" >&2
    exit 1
fi

# --- Resolve GSL_ROOT --------------------------------------------------------
: "${GSL_ROOT:=${GSL_DIR:-${GSL_HOME:-}}}"
if [[ -z "${GSL_ROOT}" ]] && command -v gsl-config >/dev/null 2>&1; then
    GSL_ROOT=$(gsl-config --prefix)
fi
if [[ -z "${GSL_ROOT}" || ! -d "${GSL_ROOT}/include" ]]; then
    echo "error: cannot locate GSL; set GSL_ROOT explicitly" >&2
    echo "       (e.g. build a local copy: ./configure --prefix=\$HOME/mysoft/gsl && make install)" >&2
    exit 1
fi

export SYSTYPE=h100 HDF5_ROOT GSL_ROOT

config=${CONFIG:-impact.conf}
if [[ ! -f "${config}" ]]; then
    echo "error: config file '${config}' not found" >&2
    exit 1
fi
cp -f "${config}" Config.sh

jobs=${JOBS:-$(nproc 2>/dev/null || echo 16)}

echo "=== planetg build ==="
echo "  config    : ${config} -> Config.sh"
echo "  systype   : h100"
echo "  HDF5_ROOT : ${HDF5_ROOT}"
echo "  GSL_ROOT  : ${GSL_ROOT}"
echo "  jobs      : ${jobs}"

make clean >/dev/null 2>&1 || true
make -j"${jobs}"

echo "=== build finished: $(ls -la GIZMO) ==="
