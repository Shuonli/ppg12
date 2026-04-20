#!/usr/bin/env bash

source /sphenix/u/shuhang98/setup.sh

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

# Merge MC. All-range configs (var_type ending in "_all") get the
# lumi-weighted cross-period merge instead of MergeSim.C, which expects
# per-sample inputs that don't exist for all-range var_types.
case "$CONFIGNAME" in
  *_all.yaml)
    PERIOD0=${CONFIGNAME/_all.yaml/_0rad.yaml}
    PERIOD1=${CONFIGNAME/_all.yaml/_1p5mrad.yaml}
    # Bare nominal "all" pairs with the historical "_nom" 1.5 mrad config.
    [[ "$CONFIGNAME" == "config_bdt_all.yaml" ]] && PERIOD1="config_bdt_nom.yaml"
    bash merge_periods.sh "$PERIOD0" "$PERIOD1" "$CONFIGNAME"
    ;;
  *)
    root -l -b -q 'MergeSim.C("'"${CONFIGNAME}"'")'
    ;;
esac

# Calculate photon yield
root -l -b -q 'CalculatePhotonYield.C("'"${CONFIGNAME}"'", true)' #isMC = true
root -l -b -q 'CalculatePhotonYield.C("'"${CONFIGNAME}"'", false)' #isMC = false
