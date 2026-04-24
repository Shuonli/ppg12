#!/usr/bin/env bash

source /sphenix/u/shuhang98/setup.sh

# Check if a config file name was provided
if [ $# -lt 1 ]; then
  echo "Usage: $0 <config_file>"
  exit 1
fi

CONFIGNAME=$1

# TTree stage: single-pass RecoEffCalculator_TTreeReader per sample + data.
# Vertex-scan Pass 1 is obsolete — all feeder configs have
# truth_vertex_reweight_on=1 which forces vertex_reweight_on=0 at Pass 2 time,
# so no *_vtxscan.root input is needed downstream.
#
# NOTE: jet5 is excluded for consistency with the DI pipeline (no jet5_double
# exists); MergeSim.C drops jet5 from the inclusive-jet merge.
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon5")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon10")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "photon20")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet8")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet12")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet20")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet30")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "jet40")' &
root -l -b -q 'RecoEffCalculator_TTreeReader.C("'"${CONFIGNAME}"'", "data")'

wait
