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

"""


import numpy as np
from .core import TrefftzMesh
from .locators import KDTreeLocator

try:
    from meshio import Mesh
except ImportError as e:
    raise ImportError(
        "The module trefftz.mesh.from_pygmsh requires pygmsh.\n"
        "Install it with: pip install trefftz[pygmsh]"
    ) from e


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

