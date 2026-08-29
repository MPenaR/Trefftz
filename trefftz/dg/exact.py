from dataclasses import dataclass
from trefftz.numpy_types import float_array, complex_array
import numpy as np

@dataclass
class PlaneWave():
    d: float_array
    A: complex
    k: float
    def __call__(self, X: float, Y: float) -> complex_array:
        return self.A*np.exp(1j*self.k*(self.d[0]*X + self.d[1]*Y))