#!/usr/bin/env bash
# Local (non-Condor) driver for variant BDT training.
# Usage: ./run_variant_training_local.sh [split|nosplit]   (default: split)
set -e
source /sphenix/u/shuhang98/setup.sh
WORKDIR="/sphenix/user/shuhangli/ppg12/FunWithxgboost"
cd "$WORKDIR"

VARIANT="${1:-split}"
FILELIST="variant_configs/filelist_${VARIANT}.txt"

if [[ ! -f "$FILELIST" ]]; then
    echo "ERROR: $FILELIST not found. Run 'python make_variant_configs.py --variant $VARIANT' first."
    exit 1
fi

while IFS= read -r cfg; do
    [[ -z "$cfg" ]] && continue
    echo "=========================================="
    echo "Training: variant_configs/$cfg"
    echo "=========================================="
    python train_variant_model.py "variant_configs/$cfg"
done < "$FILELIST"

echo "All done."
