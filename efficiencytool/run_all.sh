#!/usr/bin/env bash

# Check if a config file name was provided
if [ $# -lt 1 ]; then
  echo "Usage: $0 <config_file>"
  exit 1
fi

CONFIGNAME=$1

# Run RecoEffCalculator for the various file types
root -l -b -q 'RecoEffCalculator.C("'"${CONFIGNAME}"'", "photon5")'
root -l -b -q 'RecoEffCalculator.C("'"${CONFIGNAME}"'", "photon10")'
root -l -b -q 'RecoEffCalculator.C("'"${CONFIGNAME}"'", "photon20")'
root -l -b -q 'RecoEffCalculator.C("'"${CONFIGNAME}"'", "data")'

# Merge files
root -l -b -q 'MergeSim.C("'"${CONFIGNAME}"'")'

# Calculate photon yield
root -l -b -q 'CalculatePhotonYield.C("'"${CONFIGNAME}"'")'
