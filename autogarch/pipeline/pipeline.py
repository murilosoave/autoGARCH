import numpy as np

from typing import Union

from autogarch.estimators.mean.mean_estimator_interface import MeanEstimatorInterface
from autogarch.estimators.variance.variance_estimator_interface import VarianceEstimatorInterface
from autogarch.metrics.metric import Metric


class Pipeline:
    def __init__(self):
        self.mean_estimators: dict = {}
        self.variance_estimators: dict = {}
        self.mean_estimator: MeanEstimatorInterface = None
        self.variance_estimator: VarianceEstimatorInterface = None
        self.metrics: list = []
        self.fitted: bool = False

    def add(
        self,
        object: Union[MeanEstimatorInterface, VarianceEstimatorInterface, Metric],
        **kwargs
    ) -> "Pipeline":
        if object.type_ == "mean_estimator":
            self.mean_estimators[object] = kwargs
        elif object.type_ == "variance_estimator":
            self.variance_estimators[object] = kwargs
        elif object.type_ == "metric":
            self.metrics.append(object(**kwargs))
        else:
            raise ValueError("Invalid estimator type.")
        
        return self

    def fit(self, X: np.ndarray) -> "Pipeline":
        if not self.mean_estimators or not self.variance_estimators:
            raise RuntimeError("Pipeline must contain at least one mean and one variance estimator.")
        
        X_array = np.asarray(X, dtype=float)
        results = {"mean": {}, "variance": {}}
        for mean_estimator_class, mean_estimator_params in self.mean_estimators.items():
            for variance_estimator_class, variance_estimator_params in self.variance_estimators.items():
                mean_estimator_class.set_params(**mean_estimator_params)
                variance_estimator_class.set_params(**variance_estimator_params)
                mean_estimator_class.fit(X_array)
                variance_estimator_class.fit(mean_estimator_class.residuals)

                results["mean"][str(mean_estimator_params)] = {
                    metric_class.__class__.__name__: metric_class.calculate(mean_estimator_class.residuals)
                    for metric_class in self.metrics
                }
                results["variance"][str(variance_estimator_params)] = {
                    metric_class.__class__.__name__: metric_class.calculate(variance_estimator_class.residuals)
                    for metric_class in self.metrics
                }

                self.best_mean_estimator = min(results["mean"], key=lambda x: results["mean"][x])
                self.best_variance_estimator = min(results["variance"], key=lambda x: results["variance"][x])

                self.mean_estimator = mean_estimator_class.set_params(**mean_estimator_params)
                self.variance_estimator = variance_estimator_class.set_params(**variance_estimator_params)

        self.fitted = True

        return self
    
    def predict(self, steps: int) -> np.ndarray:
        if not self.fitted:
            raise RuntimeError("Pipeline must be fitted before making predictions.")
        
        mean_forecast = self.mean_estimator.predict(steps)
        variance_forecast = self.variance_estimator.predict(steps)

        return {"mean": mean_forecast, "variance": variance_forecast}
