from abc import ABC, abstractmethod

class Metric(ABC):
    type_: str = "metric"

    @abstractmethod
    def calculate(self, residuals):
        raise NotImplementedError
