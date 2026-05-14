"""
Utilities for constructing :class:`TrefftzMesh` objects from
``meshio`` meshes.

This module provides interoperability between ``meshio`` mesh objects
and the internal mesh representation used throughout the Trefftz library.

The main features include:

- conversion from ``meshio.Mesh`` to :class:`TrefftzMesh`
- edge reconstruction from triangle connectivity
- edge-to-triangle adjacency computation
- point location using a KD-tree accelerated search
- propagation of physical region markers from mesh files

Classes
-------
KDTreeLocator
    Spatial point locator based on triangle centroid KD-tree queries.

Functions
---------
Mesh_from_meshio
    Construct a :class:`TrefftzMesh` from a ``meshio.Mesh`` object.

Notes
-----
This module is primarily intended for meshes generated through:

- Gmsh
- pygmsh
- meshio-compatible formats

The KD-tree locator accelerates point-to-cell queries by first selecting
candidate triangles through centroid proximity and then performing
exact geometric containment checks.
"""


import numpy as np
from .core import CellLocator, TrefftzMesh
from scipy.spatial import cKDTree
from trefftz.numpy_types import float_array, int_array
from trefftz.mesh.geometry import in_triangle

try:
    from meshio import Mesh
except ImportError as e:
    raise ImportError(
        "The module trefftz.mesh.from_pygmsh requires pygmsh.\n"
        "Install it with: pip install trefftz[pygmsh]"
    ) from e


class KDTreeLocator(CellLocator):
    """
    KD-tree based locator for triangular meshes.

    This locator accelerates point-to-cell searches by building a
    KD-tree over triangle centroids and querying nearby candidate
    elements before performing exact geometric containment tests.

    Parameters
    ----------
    points : float_array
        Array of mesh vertex coordinates with shape ``(n_points, 2)``.

    triangles : int_array
        Triangle connectivity array with shape ``(n_triangles, 3)``.

    Attributes
    ----------
    tree : scipy.spatial.cKDTree
        KD-tree constructed from triangle centroids.

    radius : float
        Search radius used to retrieve candidate triangles.

    Notes
    -----
    Candidate triangles are selected based on centroid proximity and
    validated using exact point-in-triangle tests.
    """
    def __init__(self, points: float_array, triangles: int_array):
        """
        Initialize the locator and build the spatial index.
        """
        self.points = points
        self.triangles = triangles
        self.build_index()

    def build_index(self):
        """
        Construct the KD-tree spatial index.

        The tree is built from triangle centroids and an associated
        search radius is computed to guarantee retrieval of candidate
        triangles for point-location queries.
        """
        centroids = self.points[self.triangles].mean(axis=1)
        self.tree = cKDTree(centroids)
        self.radius = np.max(np.linalg.norm(self.points[self.triangles] - centroids[:, np.newaxis, :], axis=-1))

    def find_cell(self, p: float_array) -> int_array | int:
        """
        Find the triangle containing one or more query points.

        Parameters
        ----------
        p : float_array
            Query point(s).

            Accepted shapes are:

            - ``(2,)`` for a single point
            - ``(M, 2)`` for multiple points

        Returns
        -------
        int or int_array
            Index (or indices) of the containing triangle(s).

            Points outside the mesh return ``-1``.

        Raises
        ------
        ValueError
            If the input array does not have shape ``(2,)`` or ``(M, 2)``.

        Notes
        -----
        The search is performed in two stages:

        1. KD-tree retrieval of nearby candidate triangles
        2. Exact geometric containment test using
           :func:`trefftz.mesh.geometry.in_triangle`
        """
        p = np.asarray(p)
        candidates = self.tree.query_ball_point(p, self.radius)

        if p.shape == (2,):
            for i in candidates:
                if in_triangle(p, *self.points[self.triangles[i]]):
                    return i
            return -1

        elif p.ndim == 2 and p.shape[1] == 2:
            indexes = np.full(p.shape[0], dtype=np.int64, fill_value=-1)
            for j, (p_, candidates_) in enumerate(zip(p, candidates)):
                for i in candidates_:
                    if in_triangle(p_, *self.points[self.triangles[i]]):
                        indexes[j] = i
            return indexes
        else:
            raise ValueError("Input must have shape (2,) or (M, 2)")


def Mesh_from_meshio(mesh: Mesh) -> TrefftzMesh:
    """
    Construct a :class:`TrefftzMesh` from a ``meshio.Mesh`` object.

    Parameters
    ----------
    mesh : meshio.Mesh
        Input mesh containing points, triangular cells, and optional
        physical region markers.

    Returns
    -------
    TrefftzMesh
        Mesh object initialized from the meshio data structures.

    Notes
    -----
    This function performs the following operations:

    - extracts 2D point coordinates
    - reconstructs unique mesh edges from triangle connectivity
    - computes edge-to-triangle adjacency relations
    - maps boundary entities to generated edge indices
    - initializes a KD-tree based spatial locator

    Physical region markers stored in ``mesh.cell_sets_dict`` are
    propagated to the generated mesh representation.
    """
    points = mesh.points[:, 0:2]
    meshed_edges = np.sort(mesh.cells_dict["line"], axis=1)
    triangles = mesh.cells_dict["triangle"]

    # creating edges from adyacency
    edges = np.vstack([triangles[:, [0, 1]],
                       triangles[:, [1, 2]],
                       triangles[:, [2, 0]]])

    edges = np.sort(edges, axis=1)
    edges, counts = np.unique(edges, axis=0, return_counts=True)

    # pythonic loop easy to understand code, later it can be vectorized
    edge_to_index = {(i, j): idx for idx, (i, j) in enumerate(edges)}
    meshed_to_generated = np.array([edge_to_index[tuple(e)] for e in meshed_edges])

    locator = KDTreeLocator(points=points, triangles=triangles)
    cell_sets = mesh.cell_sets_dict

    for phys_ID in cell_sets.keys():
        for key in cell_sets[phys_ID].keys():
            if key == "line":
                cell_sets[phys_ID][key] = meshed_to_generated[cell_sets[phys_ID][key]]

    tri_edges = np.sort(np.stack([triangles[:, [0, 1]],
                                  triangles[:, [1, 2]],
                                  triangles[:, [2, 0]]], axis=1), axis=2)  # (T, 3, 2)

    flat_edges = tri_edges.reshape(-1, 2)       # (3T, 2)
    tri_ids = np.repeat(np.arange(len(triangles)), 3)

    # integer hashing
    max_node = len(points) # edges.max() + 1
    edge_keys = edges[:, 0].astype(np.int64) * max_node + edges[:, 1]

    # sort the keys by integer hashing in a new variable
    order = np.argsort(edge_keys)
    edge_keys_sorted = edge_keys[order]

    # flat_edges hashed
    flat_keys = flat_edges[:, 0].astype(np.int64) * max_node + flat_edges[:, 1]

    # now we can fastly search
    pos = np.searchsorted(edge_keys_sorted, flat_keys)  # position of triangle edge into the sorted global edges
    edge_ids = order[pos]  # position of the triangle edge into the global edges

    # now I want the other relation, edge to (tri1,tri2) or (tri1,-1)
    edge2triangles = np.full((len(edges), 2), -1, dtype=int)

    for E, T in zip(edge_ids, tri_ids):
        if edge2triangles[E, 0] == -1:
            edge2triangles[E, 0] = T
        else:
            edge2triangles[E, 1] = T

    return TrefftzMesh(points=points, edges=edges, triangles=triangles, edge2triangles=edge2triangles,
                locator=locator, cell_sets=cell_sets)

