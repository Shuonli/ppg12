#!/usr/bin/env bash
set -e
source /sphenix/u/shuhang98/setup.sh
WORKDIR="/sphenix/user/shuhangli/ppg12/FunWithxgboost"
cd "$WORKDIR"

while IFS= read -r cfg; do
    [[ -z "$cfg" ]] && continue
    echo "=========================================="
    echo "Training: split_configs/$cfg"
    echo "=========================================="
    python train_split_model.py "split_configs/$cfg"
done < split_configs/filelist.txt

echo "All done."
