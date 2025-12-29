from trefftz.numpy_types import float_array, int_array
import numpy as np


def in_triangle(P: float_array, A: float_array, B: float_array, C: float_array) -> bool:
    '''Computes if a point is inside a triangle'''
    AC = C - A
    AB = B - A
    AP = P - A

    u, v = np.linalg.solve(np.column_stack([AC, AB]), AP)  # computing baricentric coordinates
    tol = 1E-16

    return (u >= -tol) and (v >= -tol) and (u + v <= 1 + tol)


def triangle_area(A: float_array, B: float_array, C: float_array) -> int | int_array:
    '''Computes the area of a triangle'''
    u = (C - A).transpose()
    v = (B - A).transpose()
    det = u[0]*v[1] - u[1]*v[0]
    return 0.5*np.abs(det).transpose()
