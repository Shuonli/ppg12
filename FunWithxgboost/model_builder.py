"""
Model building utilities for binned training.
"""
from typing import Dict
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import HistGradientBoostingClassifier

try:
    import xgboost as xgb
    XGB_AVAILABLE = True
except ImportError:
    xgb = None
    XGB_AVAILABLE = False


class ModelBuilder:
    """Handles model pipeline construction and configuration."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.model_config = config['model']
        self.global_seed = config['training']['global_seed']

    def update_params(self, tuned_params: Dict) -> None:
        """Update model parameters with tuned hyperparameters.

        Args:
            tuned_params: Dictionary of hyperparameters from Optuna
        """
        if tuned_params is not None:
            self.model_config["params"].update(tuned_params)
            print(f"Updated model with tuned hyperparameters: {tuned_params}")

    def build_pipeline(self) -> Pipeline:
        """Build the complete preprocessing and model pipeline."""
        steps = []
        steps.append(("imputer", SimpleImputer(strategy="median")))
        
        if self.model_config.get("use_scaler", False):
            steps.append(("scaler", StandardScaler()))
        
        classifier = self._get_classifier()
        steps.append(("clf", classifier))
        
        return Pipeline(steps)
    
    def _get_classifier(self):
        """Get the configured classifier."""
        clf_type = self.model_config["classifier"]
        
        if clf_type == "xgb" and XGB_AVAILABLE:
            params = self.model_config["params"].copy()
            params["random_state"] = self.global_seed
            return xgb.XGBClassifier(
                **params,
                objective="binary:logistic",
                eval_metric="auc",
            )
        elif clf_type == "hgb":
            hgb_params = {
                "max_depth": None if self.model_config["params"].get("max_depth", 0) <= 0 
                           else self.model_config["params"]["max_depth"],
                "learning_rate": self.model_config["params"].get("learning_rate", 0.1),
                "max_iter": self.model_config["params"].get("n_estimators", 300),
                "l2_regularization": self.model_config["params"].get("reg_lambda", 0.0),
                "min_samples_leaf": 20,
                "random_state": self.global_seed,
            }
            return HistGradientBoostingClassifier(**hgb_params)
        else:
            raise ValueError(f"Unsupported classifier in config: {clf_type}")
    
    def get_xgb_availability(self) -> bool:
        """Check if XGBoost is available."""
        return XGB_AVAILABLE
