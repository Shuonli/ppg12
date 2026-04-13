#!/usr/bin/env bash
# Condor executable for one photon-ID BDT training job.
# Variant-agnostic: the config path is passed as $1, e.g.
#   variant_configs/config_base_v3E_split.yaml
#   variant_configs/config_base_v3E_nosplit.yaml
source /sphenix/u/shuhang98/setup.sh
WORKDIR="/sphenix/user/shuhangli/ppg12/FunWithxgboost"
cd "$WORKDIR"
echo "Training with config: $1"
python train_variant_model.py "$1"
