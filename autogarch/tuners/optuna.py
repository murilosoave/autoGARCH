import optuna

from autogarch.tuners.tuner import Tuner


class OptunaTuner(Tuner):

    def __init__(
        self,
        mean_estimator,
        variance_estimator,
        metrics,
        **kwargs,
    ):
        super().__init__(mean_estimator, variance_estimator, metrics, **kwargs)

    def tune(self, X, hiperparameter_grid, n_iters):
        def objective(trial: optuna.Trial):
            params = {
                param_name: trial.suggest_int(param_name,
                                              parameter_value[0],
                                              parameter_value[1])
                for param_name, parameter_value in hiperparameter_grid.items()
            }

            self.mean_estimator.set_params(**params)
            self.variance_estimator.set_params(**params)
            self.mean_estimator.fit(X)
            self.variance_estimator.fit(self.mean_estimator.residuals)

            metric_values = {
                metric.__class__.__name__: metric.calculate(self.mean_estimator.residuals)
                for metric in self.metrics
            }

            return metric_values["RMSE"]

        study = optuna.create_study(direction="minimize")
        study.optimize(objective, n_trials=n_iters)

        return study.best_trial.params
