#!/usr/bin/env bash

# Check if a config file name was provided
if [ $# -lt 1 ]; then
  echo "Usage: $0 <config_file>"
  exit 1
fi

CONFIGNAME=$1

# Run showershape_vertex_check for all photon (signal) samples in parallel
root -l -b -q 'showershape_vertex_check.C("'"${CONFIGNAME}"'", "photon5", true)' &
root -l -b -q 'showershape_vertex_check.C("'"${CONFIGNAME}"'", "photon10", true)' &
root -l -b -q 'showershape_vertex_check.C("'"${CONFIGNAME}"'", "photon20", true)' &
wait

cd results
hadd -f MC_vertex_check_photon.root \
  MC_efficiency_vertex_check_photon5_inclusive.root \
  MC_efficiency_vertex_check_photon10_inclusive.root \
  MC_efficiency_vertex_check_photon20_inclusive.root
