import numpy as np

from typing import Union

from autogarch.estimators.mean.mean import Mean
from autogarch.estimators.variance.variance import Variance
from autogarch.metrics.metric import Metric
from autogarch.tuners.tuner import Tuner


class Pipeline:
    def __init__(self):
        self.mean_estimators: dict = {}
        self.variance_estimators: dict = {}
        self.mean_estimator: Mean = None
        self.variance_estimator: Variance = None
        self.metrics: list = []
        self.tuner: dict[Tuner, dict] = None
        self.fitted: bool = False

    def add(
        self,
        object: Union[Mean, Variance, Metric],
        **kwargs
    ) -> "Pipeline":
        if object.type_ == "mean_estimator":
            self.mean_estimators[object] = kwargs
        elif object.type_ == "variance_estimator":
            self.variance_estimators[object] = kwargs
        elif object.type_ == "metric":
            self.metrics.append(object(**kwargs))
        elif object.type_ == "tuner":
            if not isinstance(object, Tuner):
                raise ValueError("Tuner must be an instance of Tuner.")
            if self.tuner is not None:
                raise ValueError("Pipeline can only have one tuner.")
            self.tuner = {object: kwargs}
        else:
            raise ValueError("Invalid estimator type.")
        
        return self

    def fit(self, X: np.ndarray) -> "Pipeline":
        if not self.mean_estimators or not self.variance_estimators:
            raise RuntimeError("Pipeline must contain at least one mean and one variance estimator.")
        
        X_array = np.asarray(X, dtype=float)

        tuner_class = self.tuner.popitem()[0]
        tuner_params = self.tuner.popitem()[1]

        tuner = tuner_class(**tuner_params)
        for mean_estimator_class, mean_estimator_params in self.mean_estimators.items():
            for variance_estimator_class, variance_estimator_params in self.variance_estimators.items():
                tuner.tune(
                    mean_estimator_class,
                    variance_estimator_class,
                    X_array,
                    hiperparameter_grid={
                        **mean_estimator_params,
                        **variance_estimator_params
                    },
                    n_iters=100,
                )

        return self
    
    def predict(self, steps: int) -> np.ndarray:
        if not self.fitted:
            raise RuntimeError("Pipeline must be fitted before making predictions.")
        
        mean_forecast = self.mean_estimator.predict(steps)
        variance_forecast = self.variance_estimator.predict(steps)

        return {"mean": mean_forecast, "variance": variance_forecast}
