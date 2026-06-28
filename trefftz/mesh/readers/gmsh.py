"""
Mesh readers and conversion utilities for Gmsh-generated meshes.

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


from trefftz.mesh.locators import KDTreeLocator
import numpy as np
from pathlib import Path
from trefftz.numpy_types import float_array, int_array

try:
    import gmsh
except ImportError as e:
    raise ImportError(
        "The module trefftz.mesh.from_gmsh requires gmsh.\n"
        "Install it with: pip install trefftz[gmsh]"
    ) from e


def GmshReader(file_path: Path | str) -> tuple[float_array, int_array, int_array, int_array, KDTreeLocator, dict[int, int_array]]:
    """
    Read a Gmsh(4.1.0.8) mesh file and extract mesh connectivity data.

    Parameters
    ----------
    file_path : Path or str
        Path to a Gmsh ``.msh`` file.

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

    Raises
    ------
    ValueError
        If the provided file does not have a ``.msh`` extension.

    Notes
    -----
    This function initializes the Gmsh API, loads the mesh file,
    extracts the mesh arrays using :func:`GmshArrays`, and finalizes
    the Gmsh session automatically.
    """
    file_path = Path(file_path)
    if file_path.suffix != ".msh":
        raise ValueError(f'The Gmsh Reader should be used on a GMSH generated .msh file and was called on "{file_path}" instead.')

    gmsh.initialize()
    gmsh.open(file_path.as_posix())
    points, edges, triangles, edge2triangles, locator, cell_sets = GmshArrays(model=gmsh.model)
    gmsh.finalize()
    return points, edges, triangles, edge2triangles, locator, cell_sets



def GmshArrays(model: type[gmsh.model]) -> tuple[float_array, int_array, int_array, int_array, KDTreeLocator, dict[int, int_array]]:
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
    node_tags, node_coords, node_params = model.mesh.getNodes()
    points = node_coords.reshape(-1, 3)[:, :2]  # their row-index is not valid as an ID yet
    # building look-up table  (tags not used receive an index -1 which is not valid)
    max_tag = node_tags.max()
    lut = np.full(max_tag + 1, -1)
    lut[node_tags] = np.arange(len(node_tags))


    _, tri_tags, tri_node_tags = model.mesh.getElements(dim=2)  # modify this part if the mesh becomes 3D or contains quads
    tri_node_tags = tri_node_tags[0].reshape((-1, 3))  # and this one
    triangles = lut[tri_node_tags]

    _, edge_tags, edge_node_tags = gmsh.model.mesh.getElements(dim=1)
    edge_node_tags = edge_node_tags[0].reshape(-1, 2)
    meshed_edges = np.sort(lut[edge_node_tags],axis=1)  # not sure if I need them
 
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

    phys_groups = model.getPhysicalGroups()
    phys_names = {tag: model.getPhysicalName(dim, tag) for dim, tag in phys_groups}

    phys_entities = { tag: (dim, model.getEntitiesForPhysicalGroup(dim, tag)) for dim, tag in phys_groups}
    phys_nodes = {}
    phys_tags = {}
    for tag in phys_entities:
        dim, entities = phys_entities[tag]
        phys_nodes[tag] = np.concatenate([model.mesh.getElements(dim, e)[2] for e in entities]).reshape((-1, dim+1))
        phys_tags[tag] = np.concatenate([model.mesh.getElements(dim, e)[1][0] for e in entities])-1  # watch out later for mixing quads

    cell_sets = {}
    for tag in phys_tags:
        dim, _ = phys_entities[tag]
        if dim == 1:
            cell_sets[tag] = meshed_to_generated[phys_tags[tag]]


        # for e in entities:
        #     etype, etags, enodes = model.mesh.getElements(dim, e)
        #     print(f'{enodes=}')

    
    # cell_sets = mesh.cell_sets_dict

    # for phys_ID in cell_sets.keys():
    #     for key in cell_sets[phys_ID].keys():
    #         if key == "line":
    #             cell_sets[phys_ID][key] = meshed_to_generated[cell_sets[phys_ID][key]]


    # locator = KDTreeLocator(points=points, triangles=triangles)

    return points, edges, triangles, edge2triangles, locator, cell_sets
