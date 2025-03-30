import numpy as np

from abc import ABC, abstractmethod

from autogarch.estimators.mean.mean import Mean
from autogarch.estimators.variance.variance import Variance
from autogarch.metrics.metric import Metric



class Tuner(ABC):
    type_ : str = "tuner"

    def __init__(
            self,
            mean_estimator: Mean,
            variance_estimator: Variance,
            metrics: list[Metric],
            **kwargs,
        ) -> None:
        self.mean_estimator = mean_estimator
        self.variance_estimator = variance_estimator
        self.metrics = metrics

    @abstractmethod
    def tune(self, X: np.ndarray, n_iters: int) -> dict:
        raise NotImplementedError("Tuning method not implemented.")
