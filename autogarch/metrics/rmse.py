import numpy as np

from autogarch.metrics import Metric


class RMSE(Metric):
    def __init__(self):
        super().__init__()

    def calculate(self, residuals):
        residuals_array = np.asarray(residuals, dtype=float)

        return (residuals_array ** 2).mean() ** 0.5
