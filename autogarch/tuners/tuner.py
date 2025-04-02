import numpy as np

from abc import ABC, abstractmethod

from autogarch.metrics.metric import Metric



class Tuner(ABC):
    type_ : str = "tuner"

    def __init__(
            self,
            metrics: list[Metric],
            **kwargs,
        ) -> None:
        self.metrics = metrics

    @abstractmethod
    def tune(self, X: np.ndarray, n_iters: int) -> dict:
        raise NotImplementedError("Tuning method not implemented.")
