#!/usr/bin/env bash

source /sphenix/u/shuhang98/setup.sh

# Check if a config file name was provided
if [ $# -lt 1 ]; then
  echo "Usage: $0 <config_file>"
  exit 1
fi

CONFIGNAME=$1

root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet15")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet10")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet20")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet30")' 
#root -l -b -q 'RecoEffCalculator_TTreeReader.C+("'"${CONFIGNAME}"'", "jet50")'

wait


