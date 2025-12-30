'''
Module for functions that compute geometric quantities or relations.
'''


from trefftz.numpy_types import float_array, int_array
import numpy as np


def in_triangle(P: float_array, A: float_array, B: float_array, C: float_array) -> bool:
    """
    Checks if a point is inside a given triangle.

    Parameters
    ----------
    P : (2,) float_array
        point to check.
    A : (2,) float_array
        Vertex A of triangle ABC.
    B : (2,) float_array
        Vertex B of triangle ABC.
    C : (2,) float_array
        Vertex C of triangle ABC.

    Returns
    -------
    bool
        True if the point P is inside triangle ABC.
    """
    AC = C - A
    AB = B - A
    AP = P - A

    u, v = np.linalg.solve(np.column_stack([AC, AB]), AP)  # computing baricentric coordinates
    tol = 1E-16

    return (u >= -tol) and (v >= -tol) and (u + v <= 1 + tol)


def triangle_area(A: float_array, B: float_array, C: float_array) -> int | int_array:
    """
    Computes the area of a triangle

    If several triangles are passed it computes their areas in a vectorized manner.

    Parameters
    ----------
    A : (2,) float_array
        Vertex A of triangle ABC.
    B : (2,) float_array
        Vertex B of triangle ABC.
    C : (2,) float_array
        Vertex C of triangle ABC.

    Returns
    -------
    float
        Surface area of the triangle ABC.
    """
    u = (C - A).transpose()
    v = (B - A).transpose()
    det = u[0]*v[1] - u[1]*v[0]
    return 0.5*np.abs(det).transpose()
