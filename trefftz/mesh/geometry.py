"""
Geometric utilities and spatial-query interfaces for mesh operations.

This module provides foundational geometric tools used throughout the
Trefftz library, including:

- geometric predicates
- triangle area computations
- mesh cell type definitions
- interfaces for spatial point-location algorithms

The utilities are designed to support vectorized NumPy-based workflows.

Classes
-------
CellType
    Enumeration of supported mesh cell types.

CellLocator
    Protocol defining the interface for point-to-cell location
    algorithms.

Functions
---------
in_triangle
    Test whether a point lies inside a triangle.

triangle_area
    Compute the area of one or multiple triangles.

Notes
-----
Most functions support vectorized operations over multiple geometric
entities for efficient numerical computations.
"""


from trefftz.numpy_types import float_array, int_array
import numpy as np
from typing import Protocol
from enum import IntEnum


class CellType(IntEnum):
    """
    Enumeration of supported mesh cell types.

    Attributes
    ----------
    POINT
        Zero-dimensional cell.

    SEGMENT
        One-dimensional line segment.

    TRIANGLE
        Two-dimensional triangular cell.

    Notes
    -----
    The current implementation primarily targets triangular meshes,
    although the enumeration may be extended to support additional
    element types such as quadrilaterals.
    """
    POINT = 0
    SEGMENT = 1
    TRIANGLE = 2  # consider modifying it later for quad


def in_triangle(P: float_array, A: float_array, B: float_array, C: float_array) -> bool:
    """
    Check whether a point lies inside a triangle.

    The test is performed using barycentric coordinates computed from
    the triangle vertices.

    Parameters
    ----------
    P : (2,) float_array
        Point to test.

    A : (2,) float_array
        First vertex of the triangle.

    B : (2,) float_array
        Second vertex of the triangle.

    C : (2,) float_array
        Third vertex of the triangle.

    Returns
    -------
    bool
        ``True`` if the point lies inside or on the boundary of the
        triangle, otherwise ``False``.

    Notes
    -----
    A small numerical tolerance is used to improve robustness for points
    located near triangle boundaries.
    """
    AC = C - A
    AB = B - A
    AP = P - A

    u, v = np.linalg.solve(np.column_stack([AC, AB]), AP)  # computing baricentric coordinates
    tol = 1E-16

    return (u >= -tol) and (v >= -tol) and (u + v <= 1 + tol)


def triangle_area(A: float_array, B: float_array, C: float_array) -> int | int_array:
    """
    Compute the area of one or multiple triangles.

    Parameters
    ----------
    A : float_array
        First triangle vertex or array of vertices.

    B : float_array
        Second triangle vertex or array of vertices.

    C : float_array
        Third triangle vertex or array of vertices.

    Returns
    -------
    float or float_array
        Triangle area(s).

        If vectorized inputs are provided, the returned value is an
        array containing the area of each triangle.

    Notes
    -----
    The area is computed from the determinant of the edge vectors:

    :contentReference[oaicite:0]{index=0}

    The implementation supports vectorized NumPy operations for efficient
    evaluation over multiple triangles.
    """
    u = (C - A).transpose()
    v = (B - A).transpose()
    det = u[0]*v[1] - u[1]*v[0]
    return 0.5*np.abs(det).transpose()


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
