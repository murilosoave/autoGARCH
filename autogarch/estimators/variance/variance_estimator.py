from abc import ABC, abstractmethod

from numpy.typing import ArrayLike


class VarianceEstimator(ABC):
    @abstractmethod
    def fit(self, X: ArrayLike, *args, **kwargs):
        raise NotImplementedError
    
    @abstractmethod
    def predict(self, steps: int):
        raise NotImplementedError
