#!/usr/bin/env bash

source /sphenix/u/shuhang98/setup.sh

#names=("scale20" "scale15" "scale10" "shift12" "shift15" "shift20" "bdt_base_v3E_2" "bdt_base_v3E_2" "bdt_test")
names=("bdt_none")

# Array of number pairs (low high)
number_pairs=("1 1" "2 2" "3 4" "5 8")

# Nested loops
for name in "${names[@]}"; do
    for pair in "${number_pairs[@]}"; do
        # Split the pair into low and high values
        read low high <<< "$pair"
        
        # Execute ROOT command with current parameters
        root -b -q "CONF_plots.C(\"$name\", $low, $high)" &
    done
done

