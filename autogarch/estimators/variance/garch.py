import numpy as np

from arch import arch_model
from arch.univariate import HARX
from arch.univariate.base import ARCHModelResult
from numpy.typing import ArrayLike

from .variance_estimator_interface import VarianceEstimatorInterface


class GARCH(VarianceEstimatorInterface):
    def __init__(self, p: int = 1, q: int = 1, distribution: str = "normal"):
        if p < 0 or q < 0:
            raise ValueError("Orders p and q must be non-negative integers.")
        if distribution not in ("normal", "t", "skewt"):
            raise ValueError("Distribution must be 'normal', 't' or 'skewt'.")
        self.p: int = p
        self.q: int = q
        self.distribution: str = distribution
        self.model: HARX = None
        self.model_results: ARCHModelResult = None
        self.residuals: np.ndarray = None
    
    def fit(
        self,
        X: ArrayLike,
    ) -> "GARCH":
        data = np.asarray(X, dtype=float)
        if data.ndim != 1:
            raise ValueError("Input data must be one-dimensional.")
        
        self.model = arch_model(
            data,
            p=self.p,
            q=self.q,
            dist=self.distribution
        )
        self.model_results = self.model.fit(disp="off")
        self.residuals = self.model_results.resid

        return self

    def predict(self, steps: int = 1) -> np.ndarray:
        if self.model_results is None:
            raise ValueError("Model must be fitted before forecasting.")
        if steps < 1:
            raise ValueError("Number of steps to forecast must be at least 1.")
        
        forecast = self.model_results.forecast(horizon=steps)

        return forecast.variance.values
