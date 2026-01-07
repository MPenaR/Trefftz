'''Module for defining materials, i.e. piecewise constant functions on the mesh'''

from dataclasses import dataclass
from trefftz.numpy_types import complex_array


@dataclass(slots=True, frozen=True)
class ComplexMaterial:
    values: complex_array
    mesh_id: int

    def from_regions(cls, mesh):
        pass