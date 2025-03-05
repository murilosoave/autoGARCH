from abc import ABC, abstractmethod

from numpy.typing import ArrayLike


class MeanEstimatorInterface(ABC):
    type_: str = "mean_estimator"
    input_assumptions: list = []
    output_assumptions: list = []

    @abstractmethod
    def fit(self, X: ArrayLike, *args, **kwargs):
        raise NotImplementedError
    
    @abstractmethod
    def predict(self, steps: int):
        raise NotImplementedError
    
    @abstractmethod
    def confidence_intervals(self, alpha: float):
        raise NotImplementedError
