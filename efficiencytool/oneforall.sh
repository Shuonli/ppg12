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

# Merge files
root -l -b -q 'MergeSim.C("'"${CONFIGNAME}"'")'

# Calculate photon yield
root -l -b -q 'CalculatePhotonYield.C("'"${CONFIGNAME}"'", false)' & #isMC = false
root -l -b -q 'CalculatePhotonYield.C("'"${CONFIGNAME}"'", true)' #isMC = true
