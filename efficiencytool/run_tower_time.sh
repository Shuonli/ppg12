#!/usr/bin/env bash
#set -euo pipefail

source /sphenix/u/shuhang98/setup.sh

# Run AnalyzeTruthPhotonTowers.C for photon5/photon10/photon20 and hadd outputs.
#
# Usage:
#   ./efficiencytool/run_tower_time.sh <config_file> [output_dir] [combined_output_name.root]
#
# Example:
#   ./efficiencytool/run_tower_time.sh \
#       /sphenix/user/shuhangli/ppg12/efficiencytool/config_bdt_none.yaml \
#       /sphenix/user/shuhangli/ppg12/efficiencytool/results \
#       truth_photon_tower_analysis.root

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <config_file> [output_dir] [combined_output_name.root]" >&2
  exit 1
fi

config="$1"
outdir="${2:-${repo_root}/efficiencytool/results}"
combined_name="${3:-truth_photon_tower_analysis.root}"

macro="${repo_root}/efficiencytool/AnalyzeTruthPhotonTowers.C"

if [[ ! -f "${config}" ]]; then
  echo "ERROR: config not found: ${config}" >&2
  exit 1
fi
if [[ ! -f "${macro}" ]]; then
  echo "ERROR: macro not found: ${macro}" >&2
  exit 1
fi

if ! command -v root >/dev/null 2>&1; then
  echo "ERROR: 'root' not found in PATH. Please source your ROOT environment first." >&2
  exit 1
fi
if ! command -v hadd >/dev/null 2>&1; then
  echo "ERROR: 'hadd' not found in PATH. Please source your ROOT environment first." >&2
  exit 1
fi

mkdir -p "${outdir}"

declare -a samples=(photon5 photon10 photon20)
declare -a outputs=()

echo "=================================================="
echo "Running AnalyzeTruthPhotonTowers for photon samples"
echo "=================================================="
echo "Repo root : ${repo_root}"
echo "Config    : ${config}"
echo "Out dir   : ${outdir}"
echo

for s in "${samples[@]}"; do
  outfile="${outdir}/truth_photon_tower_analysis_${s}.root"
  outputs+=("${outfile}")
  echo ">>> Running sample '${s}' -> ${outfile}"

  # Pass the whole macro call as ONE argument, like run_cluster_time.sh does.
  # Run in interpreted mode (no ACLiC compilation): note there is NO trailing '+'.
  root -l -b -q "${macro}(\"${config}\", \"${s}\", \"${outfile}\")" &
  echo
done

wait

combined="${outdir}/${combined_name}"
echo ">>> hadding into: ${combined}"
hadd -f "${combined}" "${outputs[@]}"

echo
echo "Done."
echo "  - outputs: ${outputs[*]}"
echo "  - combined: ${combined}"


