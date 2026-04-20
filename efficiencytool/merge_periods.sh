#!/usr/bin/env bash
#
# merge_periods.sh  --  run merge_periods.C on a triplet of configs
#
# Usage: merge_periods.sh <cfg_period0> <cfg_period1> <cfg_combined>
# Defaults: 0rad  + nom (1.5mrad)  ->  all
#
# Produces (paths read from cfg_combined):
#   ${eff_outfile}_${var_type}.root
#   ${eff_outfile}_jet_${var_type}.root
#   ${response_outfile}_${var_type}.root
#
# Lumi weights are derived from analysis.lumi in the per-period configs;
# the combined config's lumi field is NOT consulted for MC weighting (only
# for downstream data normalization in CalculatePhotonYield).

source /sphenix/u/shuhang98/setup.sh

CFG0=${1:-config_bdt_0rad.yaml}
CFG1=${2:-config_bdt_nom.yaml}
CFGO=${3:-config_bdt_all.yaml}

root -l -b -q "merge_periods.C(\"${CFG0}\", \"${CFG1}\", \"${CFGO}\")"
