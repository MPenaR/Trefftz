"""
Utilities for constructing and querying meshes based on
``matplotlib.tri.Triangulation``.

This module provides interoperability between Matplotlib triangulations
and the :class:`TrefftzMesh` data structures used throughout the library.

The main features include:

- conversion from Matplotlib triangulations to :class:`TrefftzMesh`
- point-to-cell location using Matplotlib's trifinder utilities
- lightweight mesh initialization for visualization and prototyping

Classes
-------
MatplotlibLocator
    Cell locator wrapper around Matplotlib triangulation search tools.

Functions
---------
Mesh_from_Matplotlib
    Construct a :class:`TrefftzMesh` from a
    ``matplotlib.tri.Triangulation``.

Notes
-----
This module is primarily intended for:

- rapid prototyping
- visualization workflows
- interoperability with Matplotlib triangulation utilities

For production mesh generation and region-aware boundary handling,
Gmsh-based workflows are generally preferred.
"""


import numpy as np
from matplotlib.tri import Triangulation
from trefftz.numpy_types import float_array, int_array
from .core import CellLocator, TrefftzMesh


class MatplotlibLocator(CellLocator):
    """
    Cell locator based on a Matplotlib triangulation.

    This class wraps the triangulation search utilities provided by
    ``matplotlib.tri.Triangulation`` and exposes a unified interface
    compatible with :class:`CellLocator`.

    Parameters
    ----------
    Tri : matplotlib.tri.Triangulation
        Triangulation object used for point-to-cell queries.

    Attributes
    ----------
    trifinder : callable
        Matplotlib triangle finder object returned by
        ``Tri.get_trifinder()``.

    Notes
    -----
    The locator uses Matplotlib's internal spatial search algorithm
    to identify the triangle containing a given point.
    """
    def __init__(self, Tri: Triangulation):
        """
        Initialize the locator from a triangulation.
        """
        self.trifinder = Tri.get_trifinder()
    
    def find_cell(self, p: float_array) -> int_array | int:
        """
        Find the triangle containing one or more points.

        Parameters
        ----------
        p : float_array
            Array of query points with shape ``(n_points, 2)`` or
            a single point of shape ``(2,)``.

        Returns
        -------
        int or int_array
            Index (or indices) of the containing triangle(s).

            Points outside the triangulation typically return ``-1``.
        """
        p_x, p_y = np.transpose(p)
        return self.trifinder(p_x, p_y)
    

def Mesh_from_Matplotlib(Tri: Triangulation) -> TrefftzMesh:
    """
    Construct a :class:`TrefftzMesh` from a Matplotlib triangulation.

    Parameters
    ----------
    Tri : matplotlib.tri.Triangulation
        Input triangulation containing points, edges, and triangle
        connectivity information.

    Returns
    -------
    TrefftzMesh
        Mesh object initialized from the triangulation data.

    Notes
    -----
    The generated mesh includes:

    - vertex coordinates
    - edge connectivity
    - triangle connectivity
    - a point locator based on Matplotlib's trifinder

    Region markers and edge-to-triangle connectivity are initialized
    as empty structures.
    """
    points = np.column_stack([Tri.x, Tri.y])
    edges = Tri.edges
    triangles = Tri.triangles
    locator = MatplotlibLocator(Tri=Tri)
    return TrefftzMesh(points=points, edges=edges, triangles=triangles, locator=locator, cell_sets={}, edge2triangles=[] )

