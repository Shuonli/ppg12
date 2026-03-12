#!/usr/bin/env bash
# hadd_isoroc.sh
# Merges per-sample roc_isoET output files into a single combined file.
# Usage: ./hadd_isoroc.sh [var_type] [results_dir]
# Example: ./hadd_isoroc.sh bdt_nom results

VAR_TYPE=${1:-bdt_nom}
OUTDIR=${2:-results}

echo "Merging var_type=${VAR_TYPE} from ${OUTDIR}/"

OUTFILE="${OUTDIR}/roc_isoET_merged_${VAR_TYPE}.root"

hadd -f "${OUTFILE}" \
  "${OUTDIR}/roc_isoET_photon5_${VAR_TYPE}.root" \
  "${OUTDIR}/roc_isoET_photon10_${VAR_TYPE}.root" \
  "${OUTDIR}/roc_isoET_photon20_${VAR_TYPE}.root" \
  "${OUTDIR}/roc_isoET_jet5_${VAR_TYPE}.root" \
  "${OUTDIR}/roc_isoET_jet12_${VAR_TYPE}.root" \
  "${OUTDIR}/roc_isoET_jet20_${VAR_TYPE}.root" \
  "${OUTDIR}/roc_isoET_jet30_${VAR_TYPE}.root" \
  "${OUTDIR}/roc_isoET_jet40_${VAR_TYPE}.root"

echo "Merged -> ${OUTFILE}"
