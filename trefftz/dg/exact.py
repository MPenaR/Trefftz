from dataclasses import dataclass
from trefftz.numpy_types import float_array


@dataclass
class PlaneWave():
    d: float_array
    A: complex