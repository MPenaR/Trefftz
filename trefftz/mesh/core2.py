"""
Core mesh data structures and geometric representations for Trefftz methods.

This module defines the primary mesh container used throughout the
Trefftz library together with the associated structured NumPy dtypes
used for efficient geometric and topological computations.

The implementation is designed around vectorized NumPy operations and
stores mesh entities as structured arrays containing precomputed
geometric quantities.

The module provides:

- triangular mesh representations
- edge and cell geometric data
- edge-to-cell connectivity
- boundary and region tagging
- normal and tangent vector construction
- spatial point-location support

Classes
-------
EdgeType
    Enumeration describing edge categories.

TrefftzMesh
    Main mesh container used by Trefftz discretizations.

Constants
---------
DIM
    Spatial dimension of the mesh representation.

edge_dtype
    Structured NumPy dtype describing edge data fields.

triangle_dtype
    Structured NumPy dtype describing triangle data fields.

Notes
-----
The mesh representation is optimized for:

- vectorized assembly routines
- discontinuous Galerkin methods
- Trefftz discretizations
- geometric queries
- boundary flux computations

Edge normals are automatically oriented consistently:

- boundary normals point outward
- inner normals follow neighboring-cell orientation conventions

The module currently targets two-dimensional triangular meshes.
"""


from typing import Final
import numpy as np
from numpy.linalg import norm
from .locators import CellLocator
from trefftz.numpy_types import float_array, int_array
from .geometry import triangle_area
from enum import Enum
from pathlib import Path


DIM: Final = 2

edge_dtype = [("P", np.float64, DIM),
              ("Q", np.float64, DIM),
              ("T", np.float64, DIM),
              ("N", np.float64, DIM),
              ("M", np.float64, DIM),
              ("l", float),
              ("on_boundary", np.bool),
              ("region", np.int8),
              ("triangles", np.int32, 2)]

arc_dtype = [("theta_1", np.float64),
                  ("theta_2", np.float64),
                  ("l", np.float64),
                  ("R", np.float64),
                  ("O", np.float64, 2),
                  ("on_boundary", np.bool),
                  ("triangles", np.int32, 2)]


def fill_edge_geometry(P: float_array, Q: float_array):
        edges = np.empty((), dtype=edge_dtype)
        edges["P"] = P
        edges["Q"] = Q
        edges["M"] = 0.5*(edges["P"]+edges["Q"])
        edges["l"] = norm(edges["Q"] - edges["P"], axis=-1)
        edges["T"] = (edges["Q"] - edges["P"])/edges["l"][..., np.newaxis]
        edges["N"] = np.column_stack([edges["T"][..., 1], -edges["T"][..., 0]])
        return edges

def fill_arc_geometry(theta_1: float, theta_2: float, R: float, center: float_array = np.array([0., 0.])):
        edges = np.empty((), dtype=arc_dtype)
        edges["R"] = R
        edges["theta_1"] = theta_1
        edges["theta_2"] = theta_2
        edges["l"]  = R*(theta_2-theta_1)
        edges["O"] = center
        return edges


triangle_dtype = [("A", np.float64, DIM),
                  ("B", np.float64, DIM),
                  ("C", np.float64, DIM),
                  ("M", np.float64, DIM),
                  ("area", np.float64)]


class TrefftzMesh[BR: Enum]:
    """
    Mesh container for Trefftz-based methods using NumPy structured arrays.

    This class stores geometric and topological information associated with
    a triangular mesh and converts it into structured NumPy arrays for
    efficient vectorized computations.

    The mesh contains:

    - Points (vertices)
    - Edges
    - Triangular cells
    - Edge-to-triangle connectivity
    - Boundary and region information

    During initialization, geometric quantities such as edge tangents,
    normals, midpoints, lengths, and triangle areas are automatically
    computed.

    Parameters
    ----------
    points : float_array
        Array of mesh vertices with shape ``(n_points, 2)``.

    edges : int_array
        Array of edge connectivity with shape ``(n_edges, 2)``.
        Each row contains the indices of the two endpoint vertices.

    triangles : int_array
        Array of triangle connectivity with shape ``(n_triangles, 3)``.
        Each row contains the indices of the three triangle vertices.

    edge2triangles : int_array
        Edge-to-triangle adjacency array with shape ``(n_edges, 2)``.
        For boundary edges, the second triangle index is ``-1``.

    locator : CellLocator
        Spatial locator used to identify the mesh cell containing a point.

    cell_sets : dict[int, int_array]
        Dictionary mapping physical region identifiers to edge indices.

    Attributes
    ----------
    edges : numpy.ndarray
        Structured array containing edge geometric and topological data.

    triangles : numpy.ndarray
        Structured array containing triangle geometric data.

    ready_for_assemble : bool
        Indicates whether all edges have a valid flux type assigned.

    Notes
    -----
    Boundary edge normals are automatically oriented outward with respect
    to the adjacent triangle.

    Inner edge normals are consistently oriented between neighboring cells.
    """

    def __init__(self, points: float_array, edges: int_array, triangles: int_array,
                 boundary_regions: type[BR], edge2triangles: int_array,
                 locator: CellLocator, cell_sets: dict[BR, int_array]):
        """
        Initialize the mesh and construct geometric data structures.
        """
        self._points = points
        self._edges = edges
        self._triangles = triangles
        self.locator = locator
        self._cell_sets = cell_sets
        self._boundary_regions = boundary_regions
        self._edge2triangles = edge2triangles
        self.construct_numpy_arrays()
        self._edges_on = cell_sets
        # I DONT LIKE THIS, MOVING TOWARDS TWO TYPES OF EDGES, and HAVE THEM SEPARATED IN THE MESH
        self.boundary_Edges = {}
        for reg in self.boundary_regions:
            self.boundary_Edges[reg] = self.edges_on(reg)


    @property
    def boundary_edges(self):
        return self.edges[self.edges["on_boundary"]]

    @property
    def interior_edges(self):
        return self.edges[~self.edges["on_boundary"]]

    def edges_on(self, region: BR):
        return self.edges[self._edges_on[region]]

    @property
    def boundary_regions(self) -> type[BR]:
        return self._boundary_regions

    def construct_numpy_arrays(self):
        """
        Construct structured NumPy arrays containing mesh geometry.

        This method computes and stores geometric quantities associated
        with edges and triangles, including:

        - Edge endpoints
        - Midpoints
        - Edge lengths
        - Tangent vectors
        - Normal vectors
        - Triangle barycenters
        - Triangle areas

        Boundary normals are oriented outward, while inner edge normals
        are consistently oriented between neighboring triangles.

        Notes
        -----
        The method updates the following attributes:

        - ``self.edges``
        - ``self.triangles``
        - ``self.ready_for_assemble``
        """
        edges = np.zeros(self.n_edges, dtype=edge_dtype)
        points = self._points
        edges["P"] = points[self._edges[:, 0], :]
        edges["Q"] = points[self._edges[:, 1], :]
        edges["M"] = 0.5*(edges["P"]+edges["Q"])
        edges["l"] = norm(edges["Q"] - edges["P"], axis=1)
        edges["T"] = 1/edges["l"][:, np.newaxis]*(edges["Q"] - edges["P"])
        edges["N"] = np.column_stack([edges["T"][:, 1], -edges["T"][:, 0]])
        edges["triangles"] = self._edge2triangles
        edges["on_boundary"] = (edges["triangles"][:, 1] == -1)
        edges["region"] = -1 # will be deprecated


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
        boundary_edges = edges[edges["on_boundary"]]
        boundary_triangles = triangles[boundary_edges["triangles"][:, 0]]
        baricenters = boundary_triangles["M"]
        midpoints = boundary_edges["M"]
        boundary_normals = np.sign(np.vecdot(midpoints-baricenters, boundary_edges["N"]))[:, np.newaxis]*boundary_edges["N"]
        edges["N"][edges["on_boundary"]] = boundary_normals

        # orienting inner normals (i don't think it should matter)

        inner_edges = edges[~edges["on_boundary"]]
        inner_triangles = triangles[inner_edges["triangles"]]
        bar_plus = inner_triangles[:, 0]["M"]
        bar_minus = inner_triangles[:, 1]["M"]

        # midpoints = boundary_edges["M"]
        inner_normals = np.sign(np.vecdot(bar_minus-bar_plus, inner_edges["N"]))[:, np.newaxis]*inner_edges["N"]
        edges["N"][~edges["on_boundary"]] = inner_normals


    def get_cell(self, p: float_array) -> int_array | int:
        """
        Find the mesh cell containing a point.

        Parameters
        ----------
        p : float_array
            Query point or points coordinates.

        Returns
        -------
        int or int_array
            Index (or indices) of the containing cell(s).
        """
        return self.locator.find_cell(p)

    @property
    def n_points(self) -> int:
        """
        Number of mesh vertices.

        Returns
        -------
        int
            Total number of points in the mesh.
        """
        return self._points.shape[0]

    @property
    def n_edges(self) -> int:
        """
        Number of mesh edges.

        Returns
        -------
        int
            Total number of edges in the mesh.
        """
        return self._edges.shape[0]
    
    @property
    def n_triangles(self) -> int:
        """
        Number of triangular cells.

        Returns
        -------
        int
            Total number of triangles in the mesh.
        """
        return self._triangles.shape[0]
    
    @classmethod
    def from_msh(cls, file_path: Path | str, boundary_regions: type[BR]) -> "TrefftzMesh[BR]":
        from trefftz.mesh.readers.gmsh import GmshReader
        """
        Create a mesh from a Gmsh .mesh file.

        Parameters
        ----------
        file_path: Path or str
            Path to the Gmsh ``.msh`` file.
        
        boundary_regions: Enum
            Enun with the names of the boundary regions

        Returns
        -------
        TrefftzMesh
            Constructed mesh instance.
        """
        points, edges, triangles, edges2triangles, locator, cell_sets = GmshReader(file_path)
        return cls(points, edges, triangles, boundary_regions, edges2triangles, locator, cell_sets)

    @classmethod
    def from_gmsh(cls, model, boundary_regions: type[BR]) -> "TrefftzMesh[BR]":
        from trefftz.mesh.readers.gmsh import GmshArrays
        """
        Create a mesh from a Gmsh model object.

        Parameters
        ----------
        model: gmsh model object (geometry has to be initialized)
            Geometry has to be initialized.

        boundary_regions: Enum
            Enun with the names of the boundary regions

        Returns
        -------
        TrefftzMesh
            Constructed mesh instance.
        """
        points, edges, triangles, edges2triangles, locator, cell_sets = GmshArrays(model)
        return cls(points, edges, triangles, boundary_regions, edges2triangles, locator, cell_sets)


    @classmethod
    def from_ngsolve(cls, mesh, boundary_regions: type[BR]) -> "TrefftzMesh[BR]":
        from trefftz.mesh.readers.ngsolve import NGsolveReader
        """
        Create a mesh from a ngsolve mesh.

        Parameters
        ----------
        mesh : NGsolve mesh
            mesh.

        boundary_regions: Enum
            Enun with the names of the boundary regions

        Returns
        -------
        TrefftzMesh
            Constructed mesh instance.
        """
        points, edges, triangles, edges2triangles, locator, cell_sets = NGsolveReader(mesh, boundary_regions)
        return cls(points, edges, triangles, boundary_regions, edges2triangles, locator, cell_sets)
    
    def curve_region(self, region: BR, radius: float, center: tuple[float, float] = (0., 0.)):
        edges = self.edges_on(region)
        ne = len(edges)
        curved_edges = np.empty(shape=(ne,), dtype=arc_dtype)
        O = np.asarray(center)
        R = radius
        for i, edge in enumerate(edges):
            P = edge["P"]
            Q = edge["Q"]

            OP = P - O
            OQ = P - Q

            theta_1 = np.atan2(OP[1], OP[0])
            theta_2 = np.atan2(OQ[1], OQ[0])
            d_theta = np.mod(theta_2 - theta_1, 2*np.pi)
            if d_theta <= np.pi:
                theta_1, theta_2 = theta_1, theta_1+d_theta
            else:
                d_theta = 2*np.pi - d_theta
                theta_1, theta_2 = theta_2, theta_2+d_theta
            
            l = R*d_theta

            curved_edges[i] = (theta_1, theta_2, l, R, O, edge["on_boundary"], edge["triangles"])

        self.boundary_Edges[region] = curved_edges