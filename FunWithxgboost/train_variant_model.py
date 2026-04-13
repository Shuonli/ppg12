import sys
from main_training import BinnedTrainingPipeline

pipeline = BinnedTrainingPipeline(sys.argv[1])
pipeline.run_training()
pipeline.evaluate_models()
pipeline.generate_plots()
pipeline.save_results()
