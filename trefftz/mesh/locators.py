"""
Spatial point-location utilities for triangular meshes.

This module defines the abstract interface for mesh cell locators
and provides a concrete implementation based on KD-tree acceleration.

Cell locators are responsible for determining which mesh element
contains a given point in space.

The module provides:

- a generic :class:`CellLocator` protocol
- a KD-tree accelerated implementation for triangular meshes
- robust point-in-triangle validation using geometric predicates

Design overview
---------------
The KD-tree approach is used to reduce the search space from the full
set of triangles to a small set of candidates based on centroid
proximity. Each candidate is then validated using an exact
point-in-triangle test.

Notes
-----
All implementations assume two-dimensional triangular meshes.
Extensions to higher dimensions or other element types may require
additional locator strategies.
"""


from trefftz.numpy_types import float_array, int_array
from typing import Protocol
from scipy.spatial import cKDTree
import numpy as np
from .geometry import in_triangle


class CellLocator(Protocol):
    """
    Protocol defining the interface for mesh cell locators.

    Cell locators are responsible for determining which mesh cell
    contains a given spatial point.

    Methods
    -------
    find_cell
        Return the index of the cell containing the query point(s).

    Notes
    -----
    Concrete implementations may rely on:

    - KD-trees
    - triangulation search structures
    - bounding volume hierarchies
    - geometric predicates
    """
    def find_cell(self, p: float_array) -> int_array | int:
        """
        Find the mesh cell containing one or more query points.

        Parameters
        ----------
        p : float_array
            Query point(s).

        Returns
        -------
        int or int_array
            Index (or indices) of the containing cell(s).
        """
        ...


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
