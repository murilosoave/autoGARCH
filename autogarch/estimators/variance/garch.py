import numpy as np

from typing import Optional

from numpy.typing import ArrayLike
from scipy.optimize import minimize

from .variance_estimator import VarianceEstimator


class GARCH(VarianceEstimator):
    def __init__(self, p: int = 1, q: int = 1):
        if p < 0 or q < 0:
            raise ValueError("Orders p and q must be non-negative integers.")
        self.p: int = p
        self.q: int = q
        self.parameters_: Optional[np.ndarray] = None
        self.conditional_volatility_: Optional[np.ndarray] = None
        self.log_likelihood_: Optional[float] = None
        # Internal storage of data (residuals) used for fitting, for forecasting use.
        self._data: Optional[np.ndarray] = None
    
    def fit(
        self,
        X: ArrayLike,
        initial_params: Optional[ArrayLike] = None, 
        tol: float = 1e-6,
        solver_options: Optional[dict] = None
    ) -> "GARCH":
        data = np.asarray(X, dtype=float)
        if data.ndim != 1:
            raise ValueError("Input data must be one-dimensional.")
        n = data.shape[0]
        # Store data internally for forecasting
        self._data = data
        # Number of parameters (omega + q alphas + p betas)
        param_count = 1 + self.q + self.p
        # Default initial guess if none provided
        if initial_params is None:
            data_var = np.var(data, ddof=1) if n > 1 else (data[0]**2 if n > 0 else 1.0)
            if data_var <= 0:
                data_var = 1e-6
            omega0 = 0.1 * data_var
            initial_guess = [omega0]
            # Initial guess for alpha coefficients
            if self.q > 0:
                alpha0 = 0.1 / self.q
                initial_guess += [alpha0] * self.q
            # Initial guess for beta coefficients
            if self.p > 0:
                beta0 = 0.8 / self.p
                initial_guess += [beta0] * self.p
            initial_params = np.array(initial_guess, dtype=float)
        else:
            initial_params = np.array(initial_params, dtype=float)
            if initial_params.size != param_count:
                raise ValueError(f"Initial parameters must have length {param_count} (got {initial_params.size}).")
        # Set bounds to enforce positivity: omega > 0, alphas >= 0, betas >= 0
        bounds = [(1e-12, None)] + [(0.0, None)] * (param_count - 1)
        # Optimize negative log-likelihood
        result = minimize(
            self._negative_log_likelihood,
            initial_params,
            args=(data,), 
            bounds=bounds,
            tol=tol,
            options=solver_options or {'disp': False}
        )
        # Store results
        self.parameters_ = result.x
        self.log_likelihood_ = -result.fun
        # Compute fitted conditional variances and volatility
        fitted_variances = self._compute_conditional_variances(data, self.parameters_)
        self.conditional_volatility_ = np.sqrt(fitted_variances)
        return self
    
    def predict(self, steps: int = 1) -> np.ndarray:
        if self.parameters_ is None or self.conditional_volatility_ is None or self._data is None:
            raise ValueError("Model must be fitted before forecasting.")
        if steps < 1:
            raise ValueError("Number of steps to forecast must be at least 1.")
        omega = self.parameters_[0]
        alphas = self.parameters_[1:1+self.q]
        betas = self.parameters_[1+self.q:]
        n = self._data.shape[0]
        # Prepare arrays for extended variances (h) and squared residuals (e2)
        h_ext = np.empty(n + steps)
        e2_ext = np.empty(n + steps)
        # Historical data
        h_ext[:n] = (self.conditional_volatility_ ** 2)
        e2_ext[:n] = self._data ** 2  # squared residuals (since mean is zero, data equals residuals)
        for t in range(n, n + steps):
            h_t = omega
            # ARCH terms (use past squared residuals: actual for t < n, expected for t >= n)
            for i in range(1, self.q + 1):
                if t - i >= 0:
                    h_t += alphas[i-1] * e2_ext[t - i]
            # GARCH terms (use past variances)
            for j in range(1, self.p + 1):
                if t - j >= 0:
                    h_t += betas[j-1] * h_ext[t - j]
            if h_t < 1e-12:
                h_t = 1e-12
            h_ext[t] = h_t
            # Expected future squared residual equals predicted variance (E[ε^2] = h)
            e2_ext[t] = h_t

        return np.sqrt(h_ext[n:n+steps])
    
    def _compute_conditional_variances(self, data: np.ndarray, params: np.ndarray) -> np.ndarray:
        n = data.shape[0]
        omega = params[0]
        alphas = params[1:1+self.q] if self.q > 0 else np.array([])
        betas = params[1+self.q:] if self.p > 0 else np.array([])
        h = np.empty(n)
        # Determine initial variance h[0]
        sum_alphas = alphas.sum() if alphas.size > 0 else 0.0
        sum_betas = betas.sum() if betas.size > 0 else 0.0
        if sum_alphas + sum_betas < 1:
            # Use unconditional variance if model is stationary
            h0 = omega / max(1.0 - (sum_alphas + sum_betas), 1e-8)
        else:
            # Use data variance as fallback for initial variance
            h0 = np.var(data, ddof=1) if n > 1 else (data[0]**2 if n > 0 else omega)
        if h0 <= 0:
            h0 = 1e-6  # ensure positivity
        h[0] = h0
        # Compute conditional variances recursively
        for t in range(1, n):
            h_t = omega
            # ARCH terms (past squared residuals)
            for i in range(1, self.q + 1):
                if t - i >= 0:
                    h_t += alphas[i-1] * (data[t - i] ** 2)
            # GARCH terms (past variances)
            for j in range(1, self.p + 1):
                if t - j >= 0:
                    h_t += betas[j-1] * h[t - j]
            # Ensure non-negativity
            if h_t < 1e-12:
                h_t = 1e-12
            h[t] = h_t
        return h
    
    def _negative_log_likelihood(self, params: np.ndarray, data: np.ndarray) -> float:
        h = self._compute_conditional_variances(data, params)
        nll = 0.5 * np.sum(np.log(2 * np.pi) + np.log(h) + (data**2) / h)

        return nll
