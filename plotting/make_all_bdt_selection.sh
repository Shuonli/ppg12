#!/usr/bin/env bash

source /sphenix/u/shuhang98/setup.sh

for config in ../efficiencytool/config_bdt*.yaml; do
  # Extract the part after '_' and before '.yaml'
  input="${config#*_}"   # removes "config_" -> gives "iso0ni0nt1t3.yaml"
  input="${input%.yaml}"  # removes the suffix ".yaml" -> gives "iso0ni0nt1t3"
  
  # Now pass the extracted input to the script
  sh make_selection_plots.sh "$input"
done
