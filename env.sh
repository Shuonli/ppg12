#! /bin/bash
# PPG12 environment. Source this file, do not execute it:
#
#     source /sphenix/user/shuhangli/ppg12/env.sh
#
# Repo-local replacement for /sphenix/u/shuhang98/setup.sh, which is not in the
# repository. It performs the same steps in the same order: sPHENIX release
# setup, MYINSTALL on the library and include paths, python venv,
# setup_local.sh, ROOT and fastjet binaries on PATH. It contains no secrets.
#
# Knobs (export before sourcing to override):
#   PPG12_SPHENIX_RELEASE  release string passed to sphenix_setup.sh   (default: new)
#   MYINSTALL              private install prefix                       (default: /sphenix/u/shuhang98/install)
#                          Must hold lib64/libyaml-cpp.so (docs/BUILD_yaml-cpp.md)
#                          and lib/libCaloAna24.so (anatreemaker, see wiki/pipeline/01-tree-making.md).
#   PPG12_VENV             python venv activated if its bin/activate exists
#                          (default: /sphenix/user/shuhangli/FMNP/conversion_venv/venv)
#
# Exports PPG12_ROOT = directory containing this file. Safe to source under
# `set -e` and `set -u`: sphenix_setup.sh reads unset variables (PGHOST at
# line 112), so nounset is suspended for the body and restored afterwards.
#
# Not migrated yet: 59 tracked macros still call
#   gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so")
# literally, so overriding MYINSTALL only takes effect for them once those
# lines read $MYINSTALL (or a bare "libyaml-cpp.so", which the LD_LIBRARY_PATH
# set here already resolves).

if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  echo "env.sh must be sourced:  source ${0}" >&2
  exit 1
fi

PPG12_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
export PPG12_ROOT

case $- in *u*) _ppg12_restore_u=1 ;; *) _ppg12_restore_u=0 ;; esac
set +u

# 1. sPHENIX software stack. "-n" first clears any previously sourced
#    sPHENIX/PHENIX variables so the release can be switched cleanly.
export PPG12_SPHENIX_RELEASE="${PPG12_SPHENIX_RELEASE:-new}"
source /opt/sphenix/core/bin/sphenix_setup.sh -n "${PPG12_SPHENIX_RELEASE}"

# 2. Private install prefix (yaml-cpp, anatreemaker library).
export MYINSTALL="${MYINSTALL:-/sphenix/u/shuhang98/install}"

export LD_LIBRARY_PATH="${MYINSTALL}/lib:${LD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="${MYINSTALL}/lib64:${LD_LIBRARY_PATH}"
export ROOT_INCLUDE_PATH="${MYINSTALL}/include:${ROOT_INCLUDE_PATH:-}"

# 3. Python venv (uproot, xgboost, ruamel.yaml, matplotlib, ...).
export PPG12_VENV="${PPG12_VENV:-/sphenix/user/shuhangli/FMNP/conversion_venv/venv}"
if [ -f "${PPG12_VENV}/bin/activate" ]; then
  # shellcheck disable=SC1091
  source "${PPG12_VENV}/bin/activate"
else
  echo "env.sh: venv ${PPG12_VENV} not found, using system python3" >&2
fi

# 4. sPHENIX helper: prepends $MYINSTALL/{lib,lib64} to LD_LIBRARY_PATH and
#    $MYINSTALL/bin to PATH, and rebuilds ROOT_INCLUDE_PATH with
#    $MYINSTALL/include and its subdirectories first (setup_root6_include_path.sh).
if [ -n "${OPT_SPHENIX:-}" ] && [ -f "${OPT_SPHENIX}/bin/setup_local.sh" ]; then
  source "${OPT_SPHENIX}/bin/setup_local.sh" "${MYINSTALL}"
else
  echo "env.sh: OPT_SPHENIX not set after sphenix_setup.sh, setup_local.sh skipped" >&2
fi

# 5. ROOT and fastjet binaries first on PATH, then per-user tools.
export FASTSYS=/opt/sphenix/core/fastjet
export PATH="${ROOTSYS}/bin:${FASTSYS}/bin:${PATH}"
if [ -d "${HOME:-/nonexistent}/.local/bin" ]; then
  export PATH="${HOME}/.local/bin:${PATH}"
fi

if [ "${_ppg12_restore_u}" = 1 ]; then set -u; fi
unset _ppg12_restore_u
return 0
