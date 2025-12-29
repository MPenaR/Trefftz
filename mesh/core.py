'''Module for importing and managing meshes'''

from typing import Final
import numpy as np
from numpy.linalg import norm
from typing import Protocol

from numpy.typing import NDArray
float_array = NDArray[np.float64]
int_array = NDArray[np.int64]


DIM: Final = 2

edge_dtype = [("P", np.float64, DIM),
              ("Q", np.float64, DIM),
              ("T", np.float64, DIM),
              ("N", np.float64, DIM),
              ("M", np.float64, DIM),
              ("l", float),
              ("boundary", bool),
              ("triangles", np.int32, 2)]


class CellLocator(Protocol):
    '''Protocol for cell locators'''
    def find_cell(self, p: float_array) -> int_array | int:
        ...


class Mesh():
    '''Holds only the relevant data
    as numpy structured-arrays for easy manipulation'''

    def __init__(self, points: float_array, edges: int_array, triangles: int_array,
                 edge2triangles: int_array,
                 locator: CellLocator, cell_sets: dict[str, dict[str, int_array]]):
        self._points = points
        self._edges = edges
        self._triangles = triangles
        self.locator = locator
        self._cell_sets = cell_sets
        self._edge2triangles = edge2triangles
        self.construct_numpy_arrays()

    def construct_numpy_arrays(self):
        edges = np.zeros(self.n_edges, dtype=edge_dtype)
        points = self._points
        edges["P"] = points[self._edges[:, 0], :]
        edges["Q"] = points[self._edges[:, 1], :]
        edges["M"] = 0.5*(edges["P"]+edges["Q"])
        edges["l"] = norm(edges["Q"] - edges["P"], axis=1)
        edges["T"] = 1/edges["l"][:, np.newaxis]*(edges["Q"] - edges["P"])
        edges["N"] = np.column_stack([edges["T"][:, 1], -edges["T"][:, 0]])
        edges["triangles"] = self._edge2triangles
        edges["boundary"] = edges["triangles"][:, 1] == -1
        self.edges = edges

    def get_cell(self, p: float_array) -> int_array | int:
        return self.locator.find_cell(p)

    @property
    def n_points(self) -> int:
        return self._points.shape[0]

    @property
    def n_edges(self) -> int:
        return self._edges.shape[0]
