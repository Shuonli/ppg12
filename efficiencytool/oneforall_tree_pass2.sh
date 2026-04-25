#!/usr/bin/env bash
# Pass-2-only re-run of RecoEffCalculator_TTreeReader for a single config.
# Assumes the _vtxscan.root files (Pass 1 output) already exist and are valid.
# Used when a transient Pass-1 data failure left an empty data vtxscan that was
# patched by copying the nominal data vtxscan file.

source /sphenix/u/shuhang98/setup.sh

if [ $# -lt 1 ]; then
  echo "Usage: $0 <config_file>"
  exit 1
fi

CONFIGNAME=$1

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
