#!/usr/bin/env bash

source /sphenix/u/shuhang98/setup.sh

for config in config_bdt*.yaml; do
  sh oneforall.sh "$config"
done
