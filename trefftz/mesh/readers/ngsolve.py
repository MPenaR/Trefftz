"""
NGsolve readers.

This module provides tools for importing meshes generated with Gmsh and
converting them into the internal :class:`TrefftzMesh` representation
used throughout the Trefftz library.

The implementation includes:

- parsing of Gmsh ``.msh`` files
- extraction of nodes, edges, and triangular cells
- construction of edge-to-triangle adjacency relations
- propagation of physical boundary tags
- KD-tree accelerated point-location utilities


Functions
---------
GmshReader
    Read a Gmsh mesh file and construct the associated mesh arrays.

GmshArrays
    Extract mesh arrays directly from an initialized Gmsh model.

Notes
-----
The current implementation targets:

- two-dimensional meshes
- triangular elements
- Gmsh version 4.x formats

Physical groups defined in Gmsh are propagated into the mesh
``cell_sets`` structure and may be used for:

- boundary condition assignment
- region tagging
- flux specification

Point-location queries are accelerated using centroid-based KD-tree
searches combined with exact geometric containment tests.
"""

from enum import StrEnum
from trefftz.mesh.locators import CellLocator
import numpy as np
# from pathlib import Path
from trefftz.numpy_types import float_array, int_array


try:
    import ngsolve
except ImportError as e:
    raise ImportError(
        "The module trefftz.mesh.from_ngsolve requires ngsolve.\n"
        "Install it with: pip install trefftz[ngsolve]"
    ) from e


class NGsolveLocator:
    def __init__(self, mesh: ngsolve.comp.Mesh):
        self._mesh = mesh
        self._elem_to_faces = np.array([E.faces[0].nr for E in mesh.Elements()])

    def find_cell(self, p: float_array) -> int_array:
        x = p[:, 0]
        y = p[:, 1]
        mp = self._mesh(x, y)
        nr = mp["nr"]
        return self._elem_to_faces[nr]


def NGsolveReader(mesh: ngsolve.comp.Mesh, boundary_regions: StrEnum) -> tuple[float_array, int_array, int_array, int_array, CellLocator, dict[int, int_array]]:
    """
    Extract mesh arrays from a Gmsh model.

    Parameters
    ----------
    model : gmsh.model
        Initialized Gmsh model containing mesh entities and physical
        groups.

    Returns
    -------
    tuple
        Tuple containing:

        - mesh vertex coordinates
        - edge connectivity
        - triangle connectivity
        - edge-to-triangle adjacency
        - spatial locator
        - physical-region cell sets

    Notes
    -----
    The extraction process includes:

    - node renumbering into compact NumPy indexing
    - edge reconstruction from triangle connectivity
    - edge uniqueness filtering
    - edge-to-triangle adjacency computation
    - physical boundary mapping
    - KD-tree locator construction

    Physical groups of dimension one are propagated into the resulting
    ``cell_sets`` dictionary and may be used for boundary-condition
    assignment.

    The current implementation assumes:

    - two-dimensional meshes
    - triangular elements
    - manifold edge connectivity
    """

    points = np.array([v.point for v in mesh.vertices])
    edges = np.array([(E.vertices[0].nr, E.vertices[1].nr) for E in mesh.edges])
    # triangles = np.array([(C.vertices[0].nr, C.vertices[1].nr, C.vertices[2].nr) for C in mesh.Elements()])
    triangles = np.array([(F.vertices[0].nr, F.vertices[1].nr, F.vertices[2].nr) for F in mesh.faces]) # face-centric? 
    # elem_to_faces = np.array([E.faces[0].nr for E in mesh.Elements()])
    edge2triangles = np.array([(E.faces[0].nr, E.faces[1].nr if len(E.faces) > 1 else -1) for E in mesh.edges])

    locator = NGsolveLocator(mesh)
    # cell_sets: dict[boundary_regions, int_array] = {}
    cell_sets = {}

    for reg in boundary_regions:
        mask = mesh.Boundaries(reg).Mask()
        cell_sets[reg] = [el.edges[0].nr for el in mesh.Elements(ngsolve.BND) if mask[el.index]]

    return points, edges, triangles, edge2triangles, locator, cell_sets
