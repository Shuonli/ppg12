#!/bin/bash
# Single-process appendix studies re-run on the current simulation trees (2026-09-10).
set -e
source /sphenix/u/shuhang98/setup.sh
cd /sphenix/user/shuhangli/ppg12/efficiencytool
case "$1" in
  deltaR) root -l -b -q 'compare_efficiency_deltaR.C("config_bdt_nom.yaml", -1)' ;;
  mbdvtx) root -l -b -q StudyMBDVertexEff.C+ ;;
  eresp)  root -l -b -q 'compare_energy_response.C("config_bdt_nom.yaml", -1)' ;;
  vtxalt) root -l -b -q VertexReweightAlt.C+ ;;
  *) echo "unknown study: $1"; exit 1 ;;
esac
