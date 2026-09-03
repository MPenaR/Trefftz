from trefftz.mesh import TrefftzMesh
from enum import Enum
from typing import Any
import numpy as np

class Elementwise[R: Enum, T]:
    def __init__(self, mesh: TrefftzMesh[Any, R],
                 parameter_mapping: dict[R, T] | None = None,
                 default: T = 1.):
        self._mesh = mesh
        self._parameter_mapping = parameter_mapping
        self._values = np.full(mesh.n_triangles,
                               fill_value=default,
                               dtype=np.asarray(list(parameter_mapping.values())).dtype)
        for reg in mesh.regions:
            self._values[mesh._reg_indexes[reg]] = parameter_mapping[reg]

    def at(self, ID: int) -> T:
        
        return np.where(ID.ravel() < 0, np.nan, self._values[ID.ravel()]).reshape(ID.shape)

    def _from_coordinates(self, x: float, y: float) -> T:
        ID = self._mesh.get_cell(x, y)
        p = self.at(ID)
        return p
    
    def __call__(self, x: float, y: float) -> T:
        return self._from_coordinates(x, y)
