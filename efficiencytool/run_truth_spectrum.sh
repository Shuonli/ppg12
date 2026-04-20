#!/usr/bin/env bash
# =============================================================================
# run_truth_spectrum.sh
# -----------------------------------------------------------------------------
# Sanity-check driver for the PPG12 cross-section pipeline. For each listed MC
# sample it runs TruthSpectrumCheck.C on combined.root and produces
# truth_spectrum_{filetype}.root in results/. At the end it runs the overlay
# macro which emits four PDFs into plotting/figures/truth_spectrum/.
#
# Samples with a missing / stale combined.root are skipped with a warning.
#
# Parallelism: processes up to PARALLEL=6 samples concurrently via xargs -P.
#
# Sample list:
#   Main pipeline signal  : photon5, photon10, photon20
#   Main pipeline BG      : jet8, jet12, jet20, jet30, jet40
#                           (see MergeSim.C:22-64; jet10/jet15/jet50 are
#                            training-input only, and jet5 has no _double
#                            partner; both are skipped here to match the merge)
#   BDT training ref      : jet5 (rendered for completeness, with lighter line
#                           style in the overlay)
#   Double-interaction    : photon{5,10,20}_double, jet{8,12,20,30,40}_double
#
# Usage:
#   bash run_truth_spectrum.sh
# =============================================================================

set -u

REPO=/sphenix/user/shuhangli/ppg12
EFF=${REPO}/efficiencytool
SIM_ROOT=${REPO}/anatreemaker/macro_maketree/sim/run28
RESULTS=${EFF}/results
FIGDIR=${REPO}/plotting/figures/truth_spectrum
PARALLEL=${PARALLEL:-6}

mkdir -p "${RESULTS}"
mkdir -p "${FIGDIR}"

SAMPLES=(
    photon5 photon10 photon20
    jet5 jet8 jet12 jet20 jet30 jet40
    photon5_double photon10_double photon20_double
    jet8_double jet12_double jet20_double jet30_double jet40_double
)

# -- Per-sample runner script emitted to a tmpfile (invoked via xargs -P).
RUNNER=$(mktemp /tmp/run_truth_one.XXXXXX.sh)
cat > "${RUNNER}" <<EOF_RUNNER
#!/bin/bash
set -u
s=\$1
EFF=${EFF}
SIM_ROOT=${SIM_ROOT}
RESULTS=${RESULTS}
input=\${SIM_ROOT}/\${s}/condorout/combined.root
if [ ! -f "\${input}" ] || [ "\$(stat -c %s "\${input}" 2>/dev/null || echo 0)" -lt 10000 ]; then
  echo "[skip \${s}] combined.root missing or stub"
  exit 0
fi
out=\${RESULTS}/truth_spectrum_\${s}.root
cd "\${EFF}" && root -l -b -q "TruthSpectrumCheck.C(\\"\${s}\\",\\"\${RESULTS}\\")" >/tmp/truth_spec_\${s}.log 2>&1
rc=\$?
if [ \$rc -ne 0 ] || [ ! -f "\${out}" ]; then
  echo "[fail \${s}] rc=\${rc}"
  exit 1
fi
echo "[ok \${s}]"
EOF_RUNNER
chmod +x "${RUNNER}"

echo "============================================================"
echo "Running TruthSpectrumCheck on ${#SAMPLES[@]} samples (parallel=${PARALLEL})"
echo "============================================================"

printf '%s\n' "${SAMPLES[@]}" | xargs -I{} -P "${PARALLEL}" bash "${RUNNER}" {}
rc=$?

rm -f "${RUNNER}"

echo
echo "============================================================"
echo "Per-sample output files"
echo "============================================================"
ls -1 "${RESULTS}"/truth_spectrum_*.root 2>/dev/null | wc -l
echo "sample results in ${RESULTS}/truth_spectrum_*.root"

echo
echo "============================================================"
echo "Running TruthSpectrumOverlay -> ${FIGDIR}"
echo "============================================================"
cd "${EFF}" && root -l -b -q "TruthSpectrumOverlay.C(\"${RESULTS}\",\"${FIGDIR}\")"

echo
echo "Done. PDFs:"
ls -1 "${FIGDIR}"/truth_spectrum_*.pdf 2>/dev/null
