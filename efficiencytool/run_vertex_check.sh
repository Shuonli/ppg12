#!/usr/bin/env bash

# Check if a config file name was provided
if [ $# -lt 1 ]; then
  echo "Usage: $0 <config_file>"
  exit 1
fi

CONFIGNAME=$1

# Run showershape_vertex_check for all jet samples in parallel
root -l -b -q 'showershape_vertex_check.C("'"${CONFIGNAME}"'", "jet10", true)' &
root -l -b -q 'showershape_vertex_check.C("'"${CONFIGNAME}"'", "jet15", true)' &
root -l -b -q 'showershape_vertex_check.C("'"${CONFIGNAME}"'", "jet20", true)' &
root -l -b -q 'showershape_vertex_check.C("'"${CONFIGNAME}"'", "jet30", true)' &
root -l -b -q 'showershape_vertex_check.C("'"${CONFIGNAME}"'", "jet50", true)' &
wait

cd results
hadd -f MC_vertex_check_jet.root \
  MC_efficiency_vertex_check_jet10_inclusive.root \
  MC_efficiency_vertex_check_jet15_inclusive.root \
  MC_efficiency_vertex_check_jet20_inclusive.root \
  MC_efficiency_vertex_check_jet30_inclusive.root \
  MC_efficiency_vertex_check_jet50_inclusive.root
