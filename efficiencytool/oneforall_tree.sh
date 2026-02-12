#!/usr/bin/env bash

source /sphenix/u/shuhang98/setup.sh

# Check if a config file name was provided
if [ $# -lt 1 ]; then
  echo "Usage: $0 <config_file>"
  exit 1
fi

CONFIGNAME=$1

# Run RecoEffCalculator for the various file types
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon5")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon10")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon20")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "data")' 

wait

