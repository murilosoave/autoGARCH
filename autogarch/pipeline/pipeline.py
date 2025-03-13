import numpy as np

from typing import Union


class Pipeline:
    def __init__(
            self,
            mean_estimator: dict[str, dict[str, Union[str, int, float]]],
            volatility_estimator: dict[str, dict[str, Union[str, int, float]]]
        ):
        self.mean_estimator = mean_estimator
        self.volatility_estimator = volatility_estimator

    def fit(self, X: np.ndarray) -> "Pipeline":
        pass
