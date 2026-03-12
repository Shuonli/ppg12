#!/usr/bin/env bash


source /sphenix/u/shuhang98/setup.sh
# Use config_showershape.yaml by default (optional first argument)
CONFIGNAME=${1:-config_showershape.yaml}
SUFFIX=${CONFIGNAME#config_}; SUFFIX=${SUFFIX%.yaml}

# Pass 1: vertex scan for all samples (fast — fills h_vertexz only)
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "photon5",  false, true)' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "photon10", false, true)' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "photon20", false, true)' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet5",     true,  true)' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet8",     true,  true)' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet12",    true,  true)' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet20",    true,  true)' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet40",    true,  true)' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "data",     false, true)'

wait

# Pass 2: full analysis for all samples
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "photon5",  false)' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "photon10", false)' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "photon20", false)' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet5",     true)' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet8",     true)' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet12",    true)' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet20",    true)' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "jet40",    true)' &
root -l -b -q 'ShowerShapeCheck.C("'"${CONFIGNAME}"'", "data",     false)'

wait


cd results
hadd -f MC_efficiencyshower_shape_signal_${SUFFIX}.root MC_efficiencyshower_shape_photon5_${SUFFIX}.root MC_efficiencyshower_shape_photon10_${SUFFIX}.root MC_efficiencyshower_shape_photon20_${SUFFIX}.root

hadd -f MC_efficiencyshower_shape_jet_inclusive_${SUFFIX}.root MC_efficiencyshower_shape_jet5_inclusive_${SUFFIX}.root MC_efficiencyshower_shape_jet8_inclusive_${SUFFIX}.root MC_efficiencyshower_shape_jet12_inclusive_${SUFFIX}.root MC_efficiencyshower_shape_jet20_inclusive_${SUFFIX}.root MC_efficiencyshower_shape_jet40_inclusive_${SUFFIX}.root
