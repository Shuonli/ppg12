#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

DATA="results/data_histo_bdt_none.root"
MC="results/MC_efficiency_bdt_none.root"
MCJET="results/MC_efficiency_jet_bdt_none.root"
OUT="results/vertex_reweight_bdt_none.root"

# Options:
#   normalize=1 (recommended for weights), rebin=0, makePdf=1
root -l -b -q "VertexReweighting.C(\"${DATA}\",\"${MC}\",\"${MCJET}\",\"${OUT}\",1,0,1)"

echo "[run_vertex_reweighting.sh] Done."
echo "  Output ROOT: ${OUT}"
echo "  Output PDF : results/vertex_reweighting.pdf"











