#!/bin/bash
# make_paper_figures.sh -- one-shot driver that re-runs every paper figure
# macro and writes outputs directly into PPG12-Paper/figures/.
#
# Usage:
#   cd plotting/paper
#   bash make_paper_figures.sh             # build all figures
#   bash make_paper_figures.sh showershape # build a single figure
#   bash make_paper_figures.sh --list      # list available targets
#
# Each target maps to one ROOT/Python macro. Default tune is bdt_nom and
# is centralized in paper_style.h (paper_tune()).

set -euo pipefail

cd "$(dirname "$0")"

PAPER_FIGDIR="../../PPG12-Paper/figures"
mkdir -p "$PAPER_FIGDIR"

# target_name -> shell command. Keep in sync with the .C / .py files in
# this directory.
declare -A TARGETS=(
    [showershape]='root -l -b -q plot_paper_showershape_yj.C'
    [iso_template]='root -l -b -q plot_paper_iso_template.C'
    [purity]='root -l -b -q plot_paper_purity.C'
    [efficiency]='root -l -b -q plot_paper_efficiency_yj.C'
    [systematics]='python plot_paper_systematics.py'
    [final]='root -l -b -q plot_paper_final_yj.C'
)

# Order matters only for log readability; macros are independent.
ORDER=(showershape iso_template purity efficiency systematics final)

run_target() {
    local name="$1"
    if [[ -z "${TARGETS[$name]+x}" ]]; then
        echo "[paper-figs] unknown target: $name" >&2
        return 2
    fi
    echo "=== [$name] ${TARGETS[$name]}"
    eval "${TARGETS[$name]}"
}

if [[ ${1:-} == "--list" ]]; then
    printf '%s\n' "${ORDER[@]}"
    exit 0
fi

if [[ $# -eq 0 ]]; then
    for t in "${ORDER[@]}"; do
        run_target "$t"
    done
else
    for t in "$@"; do
        run_target "$t"
    done
fi

echo
echo "[paper-figs] outputs written to: $(realpath "$PAPER_FIGDIR")"
ls -la "$PAPER_FIGDIR"/*.pdf 2>/dev/null | tail -15
