#!/usr/bin/env bash

source /sphenix/u/shuhang98/setup.sh

# Check if a config file name was provided
if [ $# -lt 1 ]; then
  echo "Usage: $0 <config_file>"
  exit 1
fi

CONFIGNAME=$1

# Pass 1: vertex scan for all samples (fast — fills h_vertexz only)
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon5",  true)' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon10", true)' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon20", true)' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet5",     true)' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet8",     true)' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet12",    true)' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet20",    true)' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet30",    true)' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet40",    true)' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "data",     true)'

wait

# Pass 2: full analysis for all samples
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon5")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon10")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon20")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet5")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet8")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet12")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet20")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet30")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet40")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "data")'

wait
