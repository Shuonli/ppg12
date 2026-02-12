#!/usr/bin/env bash

# Check if a config file name was provided
if [ $# -lt 1 ]; then
  echo "Usage: $0 <config_file>"
  exit 1
fi

CONFIGNAME=$1

# Run RecoEffCalculator for the various file types
#root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "photon5")'
#root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "photon10")'
#root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "photon20")'
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet10")'
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet15")'
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet30")'
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet20")'
#root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "data")'


cd results
#hadd -f MC_efficiencyshower_shape_signal.root MC_efficiencyshower_shape_photon5.root MC_efficiencyshower_shape_photon10.root MC_efficiencyshower_shape_photon20.root

hadd -f MC_efficiencyshower_shape_jet.root MC_efficiencyshower_shape_jet30.root MC_efficiencyshower_shape_jet10.root MC_efficiencyshower_shape_jet15.root MC_efficiencyshower_shape_jet20.root
