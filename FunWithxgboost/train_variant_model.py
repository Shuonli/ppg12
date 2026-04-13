"""
CLI wrapper around BinnedTrainingPipeline for a single photon-ID BDT config.

Variant-agnostic: picks up whichever config (split or nosplit) is passed on the
CLI and runs the standard training + evaluation + save pipeline. Invoked by the
Condor executable run_variant_training.sh for each entry in
variant_configs/filelist_{variant}.txt.

Usage:
    python train_variant_model.py variant_configs/config_base_v3E_split.yaml
    python train_variant_model.py variant_configs/config_base_v3E_nosplit.yaml
"""
import sys
from main_training import BinnedTrainingPipeline

pipeline = BinnedTrainingPipeline(sys.argv[1])
pipeline.run_training()
pipeline.evaluate_models()
pipeline.generate_plots()
pipeline.save_results()
