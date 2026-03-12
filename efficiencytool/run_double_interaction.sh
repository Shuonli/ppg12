#!/usr/bin/env bash

# Check if a config file name was provided
if [ $# -lt 1 ]; then
  echo "Usage: $0 <config_file>"
  exit 1
fi

CONFIGNAME=$1

# Run DoubleInteractionCheck for the various file types
root -l -b -q 'DoubleInteractionCheck.C("'"${CONFIGNAME}"'", "photon5", false)' &
root -l -b -q 'DoubleInteractionCheck.C("'"${CONFIGNAME}"'", "photon10", false)' &
root -l -b -q 'DoubleInteractionCheck.C("'"${CONFIGNAME}"'", "photon20", false)' &
root -l -b -q 'DoubleInteractionCheck.C("'"${CONFIGNAME}"'", "jet8",  true)' &
root -l -b -q 'DoubleInteractionCheck.C("'"${CONFIGNAME}"'", "jet5", true)' &
root -l -b -q 'DoubleInteractionCheck.C("'"${CONFIGNAME}"'", "jet12", true)' &
root -l -b -q 'DoubleInteractionCheck.C("'"${CONFIGNAME}"'", "jet20", true)' &
root -l -b -q 'DoubleInteractionCheck.C("'"${CONFIGNAME}"'", "jet30", true)' &
root -l -b -q 'DoubleInteractionCheck.C("'"${CONFIGNAME}"'", "jet40", true)' &
root -l -b -q 'DoubleInteractionCheck.C("'"${CONFIGNAME}"'", "data", false)'

wait

cd results
hadd -f MC_efficiency_double_interaction_check_signal.root \
  MC_efficiency_double_interaction_check_photon5.root \
  MC_efficiency_double_interaction_check_photon10.root \
  MC_efficiency_double_interaction_check_photon20.root

hadd -f MC_efficiency_double_interaction_check_jet.root \
  MC_efficiency_double_interaction_check_jet8_inclusive.root \
  MC_efficiency_double_interaction_check_jet5_inclusive.root \
  MC_efficiency_double_interaction_check_jet12_inclusive.root \
  MC_efficiency_double_interaction_check_jet20_inclusive.root \
  MC_efficiency_double_interaction_check_jet30_inclusive.root \
  MC_efficiency_double_interaction_check_jet40_inclusive.root
