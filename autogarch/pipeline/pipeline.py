import numpy as np

from typing import Union

from autogarch.estimators.mean.mean_estimator_interface import MeanEstimatorInterface
from autogarch.estimators.variance.variance_estimator_interface import VarianceEstimatorInterface


class Pipeline:
    def __init__(self):
        self.mean_estimators: list = []
        self.variance_estimators: list = []
        self.mean_estimator: MeanEstimatorInterface = None
        self.variance_estimator: VarianceEstimatorInterface = None
        self.fitted: bool = False

    def add(
        self,
        object: Union[MeanEstimatorInterface, VarianceEstimatorInterface],
        **kwargs
    ) -> "Pipeline":
        if object.type_ == "mean_estimator":
            self.mean_estimators.append({object: kwargs})
        elif object.type_ == "variance_estimator":
            self.variance_estimators.append({object: kwargs})
        else:
            raise ValueError("Invalid estimator type.")
        
        return self

    def fit(self, X: np.ndarray) -> "Pipeline":
        if not self.mean_estimators or not self.variance_estimators:
            raise RuntimeError("Pipeline must contain at least one mean and one variance estimator.")
        
        X_array = np.asarray(X, dtype=float)
        for mean_estimator in self.mean_estimators:
            for mean_estimator_class, mean_estimator_params in mean_estimator.items():
                for variance_estimator in self.variance_estimators:
                    for variance_estimator_class, variance_estimator_params in variance_estimator.items():
                        mean_estimator_class.set_params(**mean_estimator_params)
                        variance_estimator_class.set_params(**variance_estimator_params)
                        mean_estimator_class.fit(X_array)
                        variance_estimator_class.fit(mean_estimator_class.residuals)

                        self.mean_estimator = mean_estimator_class
                        self.variance_estimator = variance_estimator_class

        self.fitted = True

        return self
    
    def predict(self, steps: int) -> np.ndarray:
        if not self.fitted:
            raise RuntimeError("Pipeline must be fitted before making predictions.")
        
        mean_forecast = self.mean_estimator.predict(steps)
        variance_forecast = self.variance_estimator.predict(steps)

        return {"mean": mean_forecast, "variance": variance_forecast}
