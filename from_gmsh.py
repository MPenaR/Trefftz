import numpy as np
from .core import CellLocator, TrefftzMesh
from scipy.spatial import cKDTree
from trefftz.numpy_types import float_array, int_array
from trefftz.mesh.geometry import in_triangle

'''Module for contructing a TrefftzMesh from the gmsh API'''


try:
    import gmsh
except ImportError as e:
    raise ImportError(
        "The module trefftz.mesh.from_gmsh requires gmsh.\n"
        "Install it with: pip install trefftz[gmsh]"
    ) from e


class KDTreeLocator(CellLocator):
    def __init__(self, points: float_array, triangles: int_array):
        self.points = points
        self.triangles = triangles
        self.build_index()

    def build_index(self):
        centroids = self.points[self.triangles].mean(axis=1)
        self.tree = cKDTree(centroids)
        self.radius = np.max(np.linalg.norm(self.points[self.triangles] - centroids[:, np.newaxis, :], axis=-1))

    def find_cell(self, p: float_array) -> int_array | int:
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


def Mesh_from_gmsh(mesh: gmsh.model.mesh) -> TrefftzMesh:
    '''Returns a Mesh from a gmesh.model.mesh object'''
    node_tags, node_coords, node_params = gmsh.model.mesh.getNodes()
    points = node_coords.reshape(-1, 3)[:, :2]  # their row-index is not valid as an ID yet

    # building look-up table  (tags not used receive an index -1 which is not valid)
    max_tag = node_tags.max()
    lut = np.full(max_tag + 1, -1)
    lut[node_tags] = np.arange(len(node_tags))

# idx_connectivity = lut[tag_connectivity]

    
    
    _, tri_tags, tri_node_tags = gmsh.model.mesh.getElements(dim=2)  # modify this part if the mesh becomes 3D or contains quads
    tri_node_tags = tri_node_tags[0].reshape((-1, 3))  # and this one
    triangles = lut[tri_node_tags]

    _, edge_tags, edge_node_tags = gmsh.model.mesh.getElements(dim=1)
    edge_node_tags = edge_node_tags[0].reshape(-1, 2)
    meshed_edges = lut[edge_node_tags]  # not sure if I need them


    # meshed_edges = np.sort(mesh.cells_dict["line"], axis=1)
 
    # creating edges from adyacency
    edges = np.vstack([triangles[:, [0, 1]],
                       triangles[:, [1, 2]],
                       triangles[:, [2, 0]]])

    edges = np.sort(edges, axis=1)
    edges, counts = np.unique(edges, axis=0, return_counts=True)

    # pythonic loop easy to understand code, later it can be vectorized
    # edge_to_index = {(i, j): idx for idx, (i, j) in enumerate(edges)}
    # meshed_to_generated = np.array([edge_to_index[tuple(e)] for e in meshed_edges])

    locator = KDTreeLocator(points=points, triangles=triangles)
    # cell_sets = mesh.cell_sets_dict

    # for phys_ID in cell_sets.keys():
    #     for key in cell_sets[phys_ID].keys():
    #         if key == "line":
    #             cell_sets[phys_ID][key] = meshed_to_generated[cell_sets[phys_ID][key]]

    cell_sets = None

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
