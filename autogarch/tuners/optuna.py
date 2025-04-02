import numpy as np
import optuna

from autogarch.tuners.tuner import Tuner
from autogarch.metrics.metric import Metric
from autogarch.estimators.mean.mean import Mean
from autogarch.estimators.variance.variance import Variance


class OptunaTuner(Tuner):

    def __init__(self, metrics: list[Metric], **kwargs):
        super().__init__(metrics, **kwargs)

    def tune(
        self,
        mean_estimator: Mean,
        variance_estimator: Variance,
        X: np.ndarray,
        hiperparameter_grid: dict[str, tuple[int, int]],
        n_iters: int = 100,
    ):
        def objective(trial: optuna.Trial):
            params = {
                param_name: trial.suggest_int(param_name,
                                              parameter_value[0],
                                              parameter_value[1])
                for param_name, parameter_value in hiperparameter_grid.items()
            }

            mean_estimator.set_params(**params)
            variance_estimator.set_params(**params)

            mean_estimator.fit(X)
            variance_estimator.fit(mean_estimator.residuals)

            mean_metric_values = {
                metric.__class__.__name__: metric.calculate(mean_estimator.residuals)
                for metric in self.metrics
            }

            variance_metric_values = {
                metric.__class__.__name__: metric.calculate(variance_estimator.residuals)
                for metric in self.metrics
            }

            return (mean_metric_values["RMSE"] + variance_metric_values["RMSE"]) / 2

        study = optuna.create_study(direction="minimize")
        study.optimize(objective, n_trials=n_iters)

        return study.best_trial.params
