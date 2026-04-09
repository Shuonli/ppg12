#!/usr/bin/env bash
source /sphenix/u/shuhang98/setup.sh
WORKDIR="/sphenix/user/shuhangli/ppg12/FunWithxgboost"
cd "$WORKDIR"
echo "Training with config: $1"
python train_split_model.py "$1"
