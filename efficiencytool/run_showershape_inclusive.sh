#!/usr/bin/env bash

# Check if a config file name was provided
if [ $# -lt 1 ]; then
  echo "Usage: $0 <config_file>"
  exit 1
fi

CONFIGNAME=$1

# Run RecoEffCalculator for the various file types
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet10", true)'&
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet30", true)'&
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet20", true)'&
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet15", true)'&

wait
cd results
#hadd -f MC_efficiencyshower_shape_signal.root MC_efficiencyshower_shape_photon5.root MC_efficiencyshower_shape_photon10.root MC_efficiencyshower_shape_photon20.root

hadd -f MC_efficiencyshower_shape_jet_inclusive.root MC_efficiencyshower_shape_jet30_inclusive.root MC_efficiencyshower_shape_jet10_inclusive.root MC_efficiencyshower_shape_jet20_inclusive.root
