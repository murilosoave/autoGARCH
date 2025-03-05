import numpy as np

from math import log, pi
from typing import Union

from scipy import optimize
from scipy.special import gammaln
from scipy.stats import norm

from .mean_estimator_interface import MeanEstimatorInterface


class NormalDistribution:
    def loglikelihood(self, errors, sigma):
        n = len(errors)
        if sigma <= 0:
            return -np.inf
        log_lik = -0.5 * (n * log(2 * pi) + 2 * n * log(sigma) + np.sum((errors / sigma) ** 2))

        return log_lik


class StudentTDistribution:
    def __init__(self, df):
        self.df = df
    
    def loglikelihood(self, errors, sigma):
        nu = self.df
        n = len(errors)
        if sigma <= 0 or nu <= 0:
            return -np.inf
        term1 = gammaln((nu + 1) / 2) - gammaln(nu / 2)
        term2 = -0.5 * log(pi * nu) - log(sigma)
        sum_log = np.sum(np.log(1 + (errors / sigma) ** 2 / nu))
        log_lik = n * (term1 + term2) - (nu + 1) / 2 * sum_log

        return log_lik


class ARIMA(MeanEstimatorInterface):
    def __init__(self, p: int = 0, d: int = 0, q: int = 0, distribution: str = 'normal', df: float = None):
        self.p: int = p
        self.d: int = d
        self.q: int = q
        self.error_dist_name: str = distribution
        self.df: float = df
        if distribution == 'normal':
            self.dist: Union[NormalDistribution, StudentTDistribution] = NormalDistribution()
        elif distribution == 'student-t':
            self.dist: Union[NormalDistribution, StudentTDistribution] = StudentTDistribution(df) if df is not None else None
        else:
            raise ValueError("Unsupported error distribution. Choose 'normal' or 'student-t'.")
        self.phi = None
        self.theta = None
        self.sigma = None
        self.df_estimated = None
        self._opt_result = None

    def fit(self, X: np.ndarray) -> "ARIMA":
        series = np.asarray(X, dtype=float)
        if len(series) - self.d <= max(self.p, self.q):
            raise ValueError("Time series is too short for the specified ARIMA(p,d,q) model.")
        y_diff = self._difference_series(series)
        
        init_params = [0.1] * (self.p + self.q)
        initial_sigma = np.std(y_diff) if np.std(y_diff) > 0 else 1.0
        init_params.append(log(initial_sigma))
        if self.error_dist_name == 'student-t' and self.df is None:
            initial_nu = 10.0
            init_params.append(log(initial_nu - 2.0) if initial_nu > 2.0 else 0.0)
        init_params = np.array(init_params, dtype=float)
        
        bounds = [(-0.99, 0.99)] * (self.p + self.q) + [(None, None)]
        if self.error_dist_name == 'student-t' and self.df is None:
            bounds.append((None, None))

        result = optimize.minimize(
            self._neg_log_likelihood,
            init_params,
            args=(y_diff,),
            method='L-BFGS-B',
            bounds=bounds,
            options={"disp": True, "maxiter": 1000, "gtol": 1e-6}
        )
        
        if not result.success:
            raise RuntimeError(f"Optimization failed: {result.message}")
        
        self._opt_result = result

        opt_params = result.x
        self.phi = opt_params[:self.p] if self.p > 0 else np.array([])
        self.theta = opt_params[self.p:self.p + self.q] if self.q > 0 else np.array([])
        self.sigma = np.exp(opt_params[self.p + self.q])
        if self.error_dist_name == 'student-t':
            if self.df is None:
                raw_nu = opt_params[self.p + self.q + 1]
                self.df_estimated = 2.0 + np.exp(raw_nu)
            else:
                self.df_estimated = self.df
        else:
            self.df_estimated = None

        return self

    def predict(self, steps: int = 1) -> np.ndarray:
        if self.phi is None or self.sigma is None:
            raise RuntimeError("Model must be fitted before forecasting.")
        forecasts = []
        last_y = self._last_vals if self.d > 0 else []
        for _ in range(steps):
            ar_part = sum(self.phi[i-1] * last_y[-i] for i in range(1, self.p + 1) if len(last_y) >= i)
            ma_part = 0.0
            forecast = ar_part + ma_part
            forecasts.append(forecast)
            last_y.append(forecast)
        return np.array(forecasts)

    def confidence_intervals(self, alpha=0.05):
        if self._opt_result is None:
            raise RuntimeError("Model must be fitted before computing confidence intervals.")
        try:
            hess_inv_approx = self._opt_result.hess_inv.todense()
        except:
            hess_inv_approx = self._opt_result.hess_inv
        params_opt = self._opt_result.x
        variances = np.diag(hess_inv_approx)
        std_errors = np.sqrt(variances)

        z = norm.ppf(1 - alpha / 2)
        ci = {}

        idx = 0
        for i in range(self.p):
            est = params_opt[idx]
            se = std_errors[idx]
            ci[f"phi_{i+1}"] = (est - z*se, est + z*se)
            idx += 1
        for i in range(self.q):
            est = params_opt[idx]
            se = std_errors[idx]
            ci[f"theta_{i+1}"] = (est - z*se, est + z*se)
            idx += 1

        log_sigma_est = params_opt[idx]
        log_sigma_se = std_errors[idx]
        sigma_est = np.exp(log_sigma_est)
        sigma_se = sigma_est * log_sigma_se
        ci["sigma"] = (sigma_est - z*sigma_se, sigma_est + z*sigma_se)
        idx += 1

        if self.error_dist_name == 'student-t' and self.df is None:
            raw_nu_est = params_opt[idx]
            raw_nu_se = std_errors[idx]
            df_se = np.exp(raw_nu_est) * raw_nu_se
            ci["df"] = (
                self.df_estimated - z * df_se,
                self.df_estimated + z * df_se
            )

        return ci

    def _difference_series(self, series: np.ndarray) -> np.ndarray:
        data = np.asarray(series, dtype=float)
        last_vals = []
        if self.d > 0:
            diff_series = data
            for _ in range(self.d):
                last_vals.append(diff_series[-1])
                diff_series = np.diff(diff_series)
            y_diff = diff_series
        else:
            y_diff = data.copy()
        self._last_vals = last_vals

        return y_diff

    def _compute_residuals(self, y_diff: np.ndarray, phi: np.ndarray, theta: np.ndarray) -> np.ndarray:
        n = len(y_diff)
        residuals = np.zeros(n)
        for t in range(n):
            ar_part = sum(phi[i-1] * y_diff[t - i] for i in range(1, self.p + 1) if t - i >= 0)
            ma_part = sum(theta[j-1] * residuals[t - j] for j in range(1, self.q + 1) if t - j >= 0)
            residuals[t] = y_diff[t] - (ar_part + ma_part)

        return residuals

    def _neg_log_likelihood(self, params: list, y_diff: np.ndarray) -> float:
        p, q = self.p, self.q
        phi = params[:p] if p > 0 else np.array([])
        theta = params[p:p+q] if q > 0 else np.array([])
        log_sigma = params[p+q]
        sigma = np.exp(log_sigma)
        if self.error_dist_name == 'student-t' and self.df is None:
            raw_nu = params[p+q+1]
            nu = 2.0 + np.exp(raw_nu)
        else:
            nu = self.df
        residuals = self._compute_residuals(y_diff, phi, theta)
        if self.error_dist_name == 'normal':
            log_lik = NormalDistribution().loglikelihood(residuals, sigma)
        else:
            if self.df is None:
                n = len(residuals)
                if sigma <= 0 or nu is None or nu <= 0:
                    return np.inf
                term1 = gammaln((nu + 1) / 2) - gammaln(nu / 2)
                term2 = -0.5 * log(pi * nu) - log(sigma)
                sum_log = np.sum(np.log(1 + (residuals / sigma) ** 2 / nu))
                log_lik = n * (term1 + term2) - (nu + 1) / 2 * sum_log
            else:
                log_lik = StudentTDistribution(nu).loglikelihood(residuals, sigma)

        return -log_lik
