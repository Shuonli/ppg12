#!/usr/bin/env bash

# Use config_bdt_nom by default (optional first argument)
CONFIGNAME=${1:-config_showershape.yaml}
SUFFIX=${CONFIGNAME#config_}; SUFFIX=${SUFFIX%.yaml}

# Run RecoEffCalculator for the various file types
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "photon5")'
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "photon10")'
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "photon20")'
#root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet10")'
#root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet15")'
#root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet30")'
#root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet20")'
#root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "data")'


cd results
hadd -f MC_efficiencyshower_shape_signal_${SUFFIX}.root MC_efficiencyshower_shape_photon5_${SUFFIX}.root MC_efficiencyshower_shape_photon10_${SUFFIX}.root MC_efficiencyshower_shape_photon20_${SUFFIX}.root

#hadd -f MC_efficiencyshower_shape_jet.root MC_efficiencyshower_shape_jet30.root MC_efficiencyshower_shape_jet10.root MC_efficiencyshower_shape_jet15.root MC_efficiencyshower_shape_jet20.root
