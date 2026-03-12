#!/usr/bin/env bash

# Use config_bdt_nom by default (optional first argument)
CONFIGNAME=${1:-config_showershape.yaml}
SUFFIX=${CONFIGNAME#config_}; SUFFIX=${SUFFIX%.yaml}

# Run RecoEffCalculator for the various file types
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet5", true)'&
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet8", true)'&
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet12", true)'&
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet20", true)'&
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet30", true)'&
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet40", true)'&

wait
cd results
#hadd -f MC_efficiencyshower_shape_signal_${SUFFIX}.root MC_efficiencyshower_shape_photon5_${SUFFIX}.root MC_efficiencyshower_shape_photon10_${SUFFIX}.root MC_efficiencyshower_shape_photon20_${SUFFIX}.root

hadd -f MC_efficiencyshower_shape_jet_inclusive_${SUFFIX}.root MC_efficiencyshower_shape_jet5_inclusive_${SUFFIX}.root MC_efficiencyshower_shape_jet8_inclusive_${SUFFIX}.root MC_efficiencyshower_shape_jet12_inclusive_${SUFFIX}.root MC_efficiencyshower_shape_jet20_inclusive_${SUFFIX}.root MC_efficiencyshower_shape_jet30_inclusive_${SUFFIX}.root MC_efficiencyshower_shape_jet40_inclusive_${SUFFIX}.root
