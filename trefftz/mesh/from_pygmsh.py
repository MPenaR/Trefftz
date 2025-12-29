import numpy as np
from .core import CellLocator, Mesh
from scipy.spatial import cKDTree
from meshio import Mesh as meshioMesh
from trefftz.numpy_types import float_array, int_array


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

def Mesh_from_meshio(mesh: meshioMesh) -> Mesh:
    '''Returns a Mesh from a meshio object'''
    points = mesh.points[:, 0:2]
    meshed_edges = np.sort(mesh.cells_dict["line"], axis=1)
    triangles = mesh.cells_dict["triangle"]

    # creating edges from adyacency
    edges = np.vstack([triangles[:, [0, 1]],
                       triangles[:, [1, 2]],
                       triangles[:, [2, 0]]])

    edges = np.sort(edges, axis=1)
    edges, counts = np.unique(edges, axis=0, return_counts=True)
    boundary_edges = edges[counts == 1]
    inner_edges = edges[counts == 2]

    boundary_edges_list = np.nonzero(counts == 1)[0]
    inner_edges_list = np.nonzero(counts == 2)[0]

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

    return Mesh(points=points, edges=edges, triangles=triangles, edge2triangles=edge2triangles,
                locator=locator, cell_sets=cell_sets)

def CleanWaveGuide(R: float = 5, H: float = 1, lc: float = 0.3) -> Mesh:
    from pygmsh.geo import Geometry
    with Geometry() as geom:
        p0 = geom.add_point([-R, 0.], mesh_size=lc)
        p1 = geom.add_point([ R, 0.], mesh_size=lc)
        p2 = geom.add_point([ R, H ], mesh_size=lc)
        p3 = geom.add_point([-R, H ], mesh_size=lc)

        bottom = geom.add_line( p0, p1)
        right  = geom.add_line(p1, p2)
        top    = geom.add_line(p2, p3)
        left   = geom.add_line(p3, p0)

        boundary = geom.add_curve_loop([bottom, right, top, left])
        domain = geom.add_plane_surface(boundary)

        geom.add_physical(domain, "Omega")
        geom.add_physical([bottom, top], "Gamma")
        geom.add_physical(left, "S_L")
        geom.add_physical(right, "S_R")
        geom.add_physical([left, right], "S")

        M = geom.generate_mesh()
    return Mesh_from_meshio(M)
