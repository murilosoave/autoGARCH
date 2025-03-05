from abc import ABC, abstractmethod

from numpy.typing import ArrayLike


class VarianceEstimatorInterface(ABC):
    type_: str = "variance_estimator"
    input_assumptions: list = []
    output_assumptions: list = []

    @abstractmethod
    def fit(self, X: ArrayLike, *args, **kwargs):
        raise NotImplementedError
    
    @abstractmethod
    def predict(self, steps: int):
        raise NotImplementedError
