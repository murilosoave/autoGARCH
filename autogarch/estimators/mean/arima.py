import numpy as np

from statsmodels.tsa.arima.model import ARIMA as SM_ARIMA
from statsmodels.tsa.arima.model import ARIMAResults

from .mean import Mean


class ARIMA(Mean):
    def __init__(self, p: int = 0, d: int = 0, q: int = 0):
        self.p: int = p
        self.d: int = d
        self.q: int = q
        self.model: SM_ARIMA = None
        self.model_results: ARIMAResults = None
        self.residuals: np.ndarray = None

    def fit(self, X: np.ndarray) -> "ARIMA":
        X_array = np.asarray(X, dtype=float)
        
        self.model = SM_ARIMA(X_array, order=(self.p, self.d, self.q))
        self.model_results = self.model.fit()
        self.residuals = self.model_results.resid

        return self

    def predict(self, steps: int = 1) -> np.ndarray:
        if self.model is None:
            raise RuntimeError("Model must be fitted before making predictions.")
        forecasts = self.model_results.forecast(steps)

        return np.array(forecasts)

    def confidence_intervals(self, alpha=0.05):
        if self.model is None:
            raise RuntimeError("Model must be fitted before computing confidence intervals.")
        
        ci = {
            param_name: param_interval
            for param_name, param_interval in zip(
                self.model_results.param_names,
                self.model_results.conf_int(alpha=alpha)
            )
        }

        return ci
