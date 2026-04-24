#!/usr/bin/env bash

# Fail fast: any merge / yield step that returns nonzero aborts the script
# instead of falling through into a downstream stage that would compute on
# stale or empty data. set -u is enabled AFTER sourcing setup.sh because
# sphenix_setup.sh references unbound vars (e.g. PGHOST) which would trip -u.
set -eo pipefail

source /sphenix/u/shuhang98/setup.sh
set -u

# Check if a config file name was provided
if [ $# -lt 1 ]; then
  echo "Usage: $0 <config_file>"
  exit 1
fi

CONFIGNAME=$1

# Run RecoEffCalculator for the various file types
#oot -l -b -q 'RecoEffCalculator.C("'"${CONFIGNAME}"'", "jet15")' &
#oot -l -b -q 'RecoEffCalculator.C("'"${CONFIGNAME}"'", "jet10")' &
#oot -l -b -q 'RecoEffCalculator.C("'"${CONFIGNAME}"'", "jet20")' &
#oot -l -b -q 'RecoEffCalculator.C("'"${CONFIGNAME}"'", "jet30")' &
#oot -l -b -q 'RecoEffCalculator.C("'"${CONFIGNAME}"'", "photon5")' &
#oot -l -b -q 'RecoEffCalculator.C("'"${CONFIGNAME}"'", "photon10")' &
#oot -l -b -q 'RecoEffCalculator.C("'"${CONFIGNAME}"'", "photon20")' &
#oot -l -b -q 'RecoEffCalculator.C("'"${CONFIGNAME}"'", "data")' 

# Merge MC by dispatch on filename convention:
#   config_bdt_X_0rad.yaml    → per-period merge-feeder    → MergeSim only
#   config_bdt_X_1p5mrad.yaml → per-period merge-feeder    → MergeSim only
#   config_bdt_X.yaml (bare)  → all-range parent IF both companion merge-
#                               feeders exist on disk → merge_periods.sh +
#                               RecoEffCalculator(data) (merge_periods only
#                               handles MC, so data_histo must be generated
#                               here).
#                             → standalone (no companions) → MergeSim
#
# The per-period MergeSim for merge-feeders is redundant when Phase 1
# (oneforall_tree_double.sh) already ran it, but MergeSim is idempotent.
case "$CONFIGNAME" in
  *_0rad.yaml|*_1p5mrad.yaml)
    root -l -b -q 'MergeSim.C("'"${CONFIGNAME}"'")'
    ;;
  *)
    base=${CONFIGNAME%.yaml}
    PERIOD0="${base}_0rad.yaml"
    PERIOD1="${base}_1p5mrad.yaml"
    if [ -f "$PERIOD0" ] && [ -f "$PERIOD1" ]; then
      bash merge_periods.sh "$PERIOD0" "$PERIOD1" "$CONFIGNAME"
      # Data was processed per-period in Phase 1; hadd to build all-range.
      # Replaces the earlier RecoEffCalculator_TTreeReader(..., "data") re-read,
      # which was redundant (per-period run ranges are disjoint + their union =
      # all-range, and data histograms are unweighted, so hadd is exact).
      var_all=$(grep -E "^  var_type:" "$CONFIGNAME" | awk -F'"' '{print $2}')
      var_0=$(grep -E "^  var_type:"   "$PERIOD0"    | awk -F'"' '{print $2}')
      var_1=$(grep -E "^  var_type:"   "$PERIOD1"    | awk -F'"' '{print $2}')
      hadd -f "results/data_histo_${var_all}.root" \
              "results/data_histo_${var_0}.root" \
              "results/data_histo_${var_1}.root"
    else
      # Standalone config with no per-period siblings: legacy MergeSim path.
      root -l -b -q 'MergeSim.C("'"${CONFIGNAME}"'")'
    fi
    ;;
esac

# Calculate photon yield
root -l -b -q 'CalculatePhotonYield.C("'"${CONFIGNAME}"'", true)' #isMC = true
root -l -b -q 'CalculatePhotonYield.C("'"${CONFIGNAME}"'", false)' #isMC = false

# Selection plots (sPHENIX-style PDFs in plotting/figures/). Bare-name configs
# only — per-period feeders don't get their own plot set (plots for nominal
# periods can be generated manually via make_selection_plots.sh bdt_nom_0rad).
case "$CONFIGNAME" in
  *_0rad.yaml|*_1p5mrad.yaml|*_0mrad.yaml)
    ;;
  *)
    var_type=$(grep -E "^  var_type:" "$CONFIGNAME" | awk -F'"' '{print $2}')
    bash /sphenix/user/shuhangli/ppg12/plotting/make_selection_plots.sh "$var_type"
    ;;
esac
