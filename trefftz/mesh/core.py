'''Module for importing and managing meshes'''

from typing import Final
import numpy as np
from numpy.linalg import norm
from .geometry import CellLocator
from trefftz.numpy_types import float_array, int_array
from .geometry import triangle_area
from enum import IntEnum
from pathlib import Path
from .readers import GmshReader
from .geometry import CellType

class EdgeType(IntEnum):
    INNER = 0
    BOUNDARY = 1

DIM: Final = 2

edge_dtype = [("P", np.float64, DIM),
              ("Q", np.float64, DIM),
              ("T", np.float64, DIM),
              ("N", np.float64, DIM),
              ("M", np.float64, DIM),
              ("l", float),
              ("type", np.int8),
              ("flux_type", np.int8),
              ("region", np.int8),
              ("triangles", np.int32, 2)]

triangle_dtype = [("A", np.float64, DIM),
                  ("B", np.float64, DIM),
                  ("C", np.float64, DIM),
                  ("M", np.float64, DIM),
                  ("area", np.float64)]

class TrefftzMesh():
    '''Returns only the relevant data
    as numpy structured-arrays for easy manipulation'''

    def __init__(self, points: float_array, edges: int_array, triangles: int_array,
                 edge2triangles: int_array,
                 locator: CellLocator, cell_sets: dict[int, int_array]):
        self._points = points
        self._edges = edges
        self._triangles = triangles
        self.locator = locator
        self._cell_sets = cell_sets
        self._edge2triangles = edge2triangles
        self.ready_for_assemble = False
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
        edges["type"] = (edges["triangles"][:, 1] == -1).astype(np.int8)
        edges["flux_type"] = -1
        # edges["flux_type"][edges["type"] == EdgeType.INNER] = FluxType.TRANSMISSION
        edges["flux_type"][edges["type"] == EdgeType.INNER] = 0
        edges["region"] = -1
        cell_sets_1D = self._cell_sets
        for region in cell_sets_1D:
            edges["region"][cell_sets_1D[region]] = region

        self.edges = edges

        triangles = np.zeros(self.n_triangles, dtype=triangle_dtype)
        triangles["A"] = points[self._triangles[:, 0], :]
        triangles["B"] = points[self._triangles[:, 1], :]
        triangles["C"] = points[self._triangles[:, 2], :]
        triangles["M"] = 1/3*(triangles["A"] + triangles["B"] + triangles["C"])
        triangles["area"] = triangle_area(A=triangles["A"],
                                          B=triangles["B"],
                                          C=triangles["C"])        
        self.triangles = triangles

        # orienting boundary normals
        boundary_edges = edges[edges["type"] == EdgeType.BOUNDARY]
        boundary_triangles = triangles[boundary_edges["triangles"][:, 0]]
        baricenters = boundary_triangles["M"]
        midpoints = boundary_edges["M"]
        boundary_normals = np.sign(np.vecdot(midpoints-baricenters, boundary_edges["N"]))[:, np.newaxis]*boundary_edges["N"]
        edges["N"][edges["type"] == EdgeType.BOUNDARY] = boundary_normals

        # orienting inner normals (i don't think it should matter)

        inner_edges = edges[edges["type"] == EdgeType.INNER]
        inner_triangles = triangles[inner_edges["triangles"]]
        bar_plus = inner_triangles[:, 0]["M"]
        bar_minus = inner_triangles[:, 1]["M"]
        
        #midpoints = boundary_edges["M"]
        inner_normals = np.sign(np.vecdot(bar_minus-bar_plus, inner_edges["N"]))[:, np.newaxis]*inner_edges["N"]
        edges["N"][edges["type"] == EdgeType.INNER] = inner_normals

        if np.all(self.edges["flux_type"] >= 0):
            print('Mesh is ready to be assebled')
            self.ready_for_assemble = True

    def get_cell(self, p: float_array) -> int_array | int:
        return self.locator.find_cell(p)

    @property
    def n_points(self) -> int:
        return self._points.shape[0]

    @property
    def n_edges(self) -> int:
        return self._edges.shape[0]
    
    @property
    def n_triangles(self) -> int:
        return self._triangles.shape[0]
    
    @classmethod
    def from_gmsh(cls, file_path: Path | str):
        points, edges, triangles, edges2triangles, locator, cell_sets = GmshReader(file_path)
        return cls(points, edges, triangles, edges2triangles, locator, cell_sets)
    
            
