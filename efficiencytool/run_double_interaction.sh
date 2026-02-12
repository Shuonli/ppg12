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
root -l -b -q 'DoubleInteractionCheck.C("'"${CONFIGNAME}"'", "jet10", false)' &
root -l -b -q 'DoubleInteractionCheck.C("'"${CONFIGNAME}"'", "jet15", false)' &
root -l -b -q 'DoubleInteractionCheck.C("'"${CONFIGNAME}"'", "jet20", false)' &
root -l -b -q 'DoubleInteractionCheck.C("'"${CONFIGNAME}"'", "jet30", false)' &
root -l -b -q 'DoubleInteractionCheck.C("'"${CONFIGNAME}"'", "data", false)'

wait

cd results
hadd -f MC_efficiency_double_interaction_check_signal.root \
  MC_efficiency_double_interaction_check_photon5.root \
  MC_efficiency_double_interaction_check_photon10.root \
  MC_efficiency_double_interaction_check_photon20.root

hadd -f MC_efficiency_double_interaction_check_jet.root \
  MC_efficiency_double_interaction_check_jet10.root \
  MC_efficiency_double_interaction_check_jet15.root \
  MC_efficiency_double_interaction_check_jet20.root \
  MC_efficiency_double_interaction_check_jet30.root
