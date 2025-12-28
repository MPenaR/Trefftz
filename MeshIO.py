'''Module for importing and managing meshes'''

from typing import Final
import numpy as np
from numpy.linalg import norm
from matplotlib.tri import Triangulation
import matplotlib.pyplot as plt 
from numpy.typing import NDArray
from scipy.spatial import cKDTree
from typing import Protocol
from meshio import Mesh as meshioMesh

float_array = NDArray[np.float64]
int_array = NDArray[np.int64]

# import ngsolve as ns

DIM: Final = 2

edge_dtype = [("P", np.float64, DIM),
              ("Q", np.float64, DIM),
              ("T", np.float64, DIM),
              ("N", np.float64, DIM),
              ("M", np.float64, DIM),
              ("l", float),
              ("boundary", bool),
              ("triangles", np.int32, 2)]


def in_triangle(P: float | float_array, A: float | float_array, B: float | float_array, C: float | float_array) -> bool:
    '''Computes if a point is inside a triangle'''
    AC = C - A
    AB = B - A
    AP = P - A

    u, v = np.linalg.solve(np.column_stack([AC, AB]), AP)  # computing baricentric coordinates
    tol = 1E-16

    return (u >= -tol) and (v >= -tol) and (u + v <= 1 + tol)


class _CellLocator(Protocol):
    '''Protocol for cell locators'''
    def find_cell(self, p: float_array) -> int_array | int:
        ...


class _MatplotlibLocator(_CellLocator):
    def __init__(self, Tri: Triangulation):
        self.trifinder = Tri.get_trifinder()
    
    def find_cell(self, p: float_array) -> int_array | int:
        p_x, p_y = np.transpose(p)
        return self.trifinder(p_x, p_y)
    

class _KDTreeLocator(_CellLocator):
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


# def construct_edge_connectivity(self):
    
#     edges = np.array([frozenset(E) for E in self._edges])

#     triangles = np.array([frozenset(T) for T in self._triangles])
        
#     numpy_subset = np.frompyfunc(lambda A, B: A <= B, 2, 1)
#     mask = numpy_subset(edges[:,None], triangles[None,:])
#     self.boundary_edges_list = np.where(mask.sum(axis=1)==1)[0]
#     self.inner_edges_list = np.where(mask.sum(axis=1)==2)[0]


class Mesh():
    '''Holds only the relevant data
    as numpy structured-arrays for easy manipulation'''

    def __init__(self, points: float_array, edges: int_array, triangles: int_array,
                 edges2triangles: int_array,
                 locator: _CellLocator, cell_sets: dict[str, dict[str, int_array]]):
        self._points = points
        self._edges = edges
        self._triangles = triangles
        self.locator = locator
        self._cell_sets = cell_sets
        self._edges2triangles = edges2triangles
        self.boundary_edges_list = np.nonzero(edges2triangles[:, 1] == -1)[0]
        self.construct_numpy_arrays()

    def construct_numpy_arrays(self):
        edges = np.zeros(self.n_edges, dtype=edge_dtype)
        points = self._points
        edges["P"] = points[self._edges[:, 0], :]
        edges["Q"] = points[self._edges[:, 1], :]
        edges["M"] = 0.5*(edges["P"]+edges["Q"])
        edges["l"] = norm(edges["Q"] - edges["P"], axis=1)
        edges["T"] = 1/edges["l"][:, np.newaxis]*(edges["Q"] - edges["P"])
        edges["N"] = np.column_stack([edges["T"][:, 1], -edges["T"][:, 0]])
        edges["boundary"] = False
        edges["boundary"][self.boundary_edges_list] = True # if I change the order, edges[I] is a copy, not a view
        edges["triangles"] = self._edges2triangles
        self.edges = edges


    def get_cell(self, p: float_array) -> int_array | int:
        return self.locator.find_cell(p)

    @property
    def n_points(self) -> int:
        return self._points.shape[0]

    @property
    def n_edges(self) -> int:
        return self._edges.shape[0]
    # def plot_mesh(self):
    #     plt.triplot(Triangulation(x=self._points[:, 0], y=self._points[:, 1], triangles=self._triangles))
    #     plt.show()


def Mesh_from_Matplotlib(Tri: Triangulation):
    '''Returns a Mesh from a Matptlotlib.tri.Triangulation'''
    points = np.column_stack([Tri.x, Tri.y])
    edges = Tri.edges
    triangles = Tri.triangles
    locator = _MatplotlibLocator(Tri=Tri)
    return Mesh(points=points, edges=edges, triangles=triangles, locator=locator)


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

    locator = _KDTreeLocator(points=points, triangles=triangles)
    cell_sets = mesh.cell_sets_dict

    for phys_ID in cell_sets.keys():
        for key in cell_sets[phys_ID].keys():
            if key == "line":
                cell_sets[phys_ID][key] = meshed_to_generated[cell_sets[phys_ID][key]]


# # old method
#     tri_edges = np.sort( np.stack([triangles[:, [0, 1]],
#                                    triangles[:, [1, 2]],
#                                    triangles[:, [2, 0]],], axis=1), axis=2)  # shape (nT, 3, 2)

#     flat_edges = tri_edges.reshape(-1, 2)   # (3*nT, 2)
#     tri_ids = np.repeat(np.arange(len(triangles)), 3)




#     # Building lookup from flat_edges
#     # interior edges
#     edge_keys_inner: dict[tuple[int, int], list[int]] = {tuple(e): [] for e in inner_edges}

#     for e, t in zip(flat_edges, tri_ids):
#         key = tuple(e)
#         if key in edge_keys_inner:
#             edge_keys_inner[key].append(t)

#     inner_edges_triangles = np.array([edge_keys_inner[tuple(e)] for e in inner_edges], dtype=int)  # (nE_int, 2)

#     # boundary edges
#     edge_keys_boundary: dict[tuple[int, int], int] = {tuple(e): -1 for e in boundary_edges}

#     for e, t in zip(flat_edges, tri_ids):
#         key = tuple(e)
#         if key in edge_keys_boundary:
#             edge_keys_boundary[key] = t

#     boundary_edges_triangle = np.array([edge_keys_boundary[tuple(e)] for e in boundary_edges], dtype=int)  # (nE_bnd,)

#new method

    tri_edges = np.sort( np.stack([triangles[:, [0, 1]],
                                   triangles[:, [1, 2]],
                                   triangles[:, [2, 0]]], axis=1), axis=2)  # (T, 3, 2)

    flat_edges = tri_edges.reshape(-1, 2)       # (3T, 2)
    tri_ids = np.repeat(np.arange(len(triangles)), 3)

    # integer hashing
    max_node = len(points) # edges.max() + 1
    edge_keys = edges[:, 0] * max_node + edges[:, 1]

    # sort the keys by integer hashing in a new variable
    order = np.argsort(edge_keys)
    edge_keys_sorted = edge_keys[order]

    # flat_edges hashed
    flat_keys = flat_edges[:, 0] * max_node + flat_edges[:, 1]

    # now we can fastly search
    pos = np.searchsorted(edge_keys_sorted, flat_keys) # position of triangle edge into the sorted global edges
    edge_ids = order[pos] # position of the triangle edge into the global edges

    #now I want the other relation, edge to (tri1,tri2) or (tri1,-1)
    edge2triangles = np.full((len(edges), 2), -1, dtype=int)

    for E, T in zip(edge_ids, tri_ids):
        if edge2triangles[E, 0] == -1:
            edge2triangles[E, 0] = T
        else:
            edge2triangles[E, 1] = T


    return Mesh(points=points, edges=edges, triangles=triangles, edges2triangles=edge2triangles,
                locator=locator, cell_sets=cell_sets)


# def triangle_id_tester(M: Mesh):
#     fig, ax = plt.subplots()
#     xmin, xmax = np.min(M._points[:,0]), np.max(M._points[:,0]) 
#     ymin, ymax = np.min(M._points[:,1]), np.max(M._points[:,1])
#     N = 200
#     x = np.linspace(xmin,xmax,N)
#     y = np.linspace(ymin,ymax,N)
#     X, Y = np.meshgrid(x,y)
#     xy = np.column_stack([X.flatten(), Y.flatten()])
#     Z = mesh.get_cell(xy).reshape(X.shape).astype(float)
#     Z[Z==-1] = np.nan
#     pc = ax.pcolormesh(X,Y,Z)
#     lw = 2
#     ax.triplot(Triangulation(x=M._points[:,0], y=M._points[:,1], triangles=M._triangles),linewidth=lw, color='k')
#     for e_ID in M.boundary_edges_list:
#         P, Q = M._edges[e_ID]
#         p_x, p_y = M._points[P, :]
#         q_x, q_y = M._points[Q, :]
#         ax.plot([p_x,q_x],[p_y,q_y],'r',linewidth=lw)
#     fig.colorbar(mappable=pc)
#     ax.axis('equal')
#     def hover_text(x, y):
#         return f'ID = {M.get_cell([x, y])}'
#     annot = ax.annotate(
#         "",
#         xy=(0, 0),
#         xytext=(10, 10),
#         textcoords="offset points",
#         bbox=dict(boxstyle="round", fc="w"),
#         arrowprops=dict(arrowstyle="->"),
#     )
#     annot.set_visible(False)
#     annot.arrowprops = None
#     def on_move(event):
#         if event.inaxes != ax:
#             annot.set_visible(False)
#             fig.canvas.draw_idle()
#             return
#         xdata, ydata = event.xdata, event.ydata
#         # Update annotation
#         annot.xy = (xdata, ydata)
#         annot.set_text(hover_text(xdata, ydata))
#         annot.set_visible(True)
#         fig.canvas.draw_idle()

#     # Connect event
#     fig.canvas.mpl_connect("motion_notify_event", on_move)
#     plt.show()


def plot_waveguide(M: Mesh, plot_tangents: bool = False, plot_normals: bool = False):
    from matplotlib.collections import LineCollection
    _, ax = plt.subplots()

    lw = 1
    # ax.triplot(Triangulation(x=M._points[:,0], y=M._points[:,1], triangles=M._triangles),linewidth=lw, color='k')

    S = M._cell_sets["S"]["line"]
    G = M._cell_sets["Gamma"]["line"]  # it allows for multidimensional subsets

    inner = np.where(np.logical_not(M.edges["boundary"]))[0]

    ax.add_collection(LineCollection(np.stack([mesh.edges[inner]["P"], mesh.edges[inner]["Q"]], axis=1), 
                                     colors='k', linewidths=lw))
    ax.add_collection(LineCollection(np.stack([mesh.edges[S]["P"], mesh.edges[S]["Q"]], axis=1), 
                                     colors='r', linewidths=lw))
    ax.add_collection(LineCollection(np.stack([mesh.edges[G]["P"], mesh.edges[G]["Q"]], axis=1), 
                                     colors='b', linewidths=lw))
    
    if plot_tangents:
        ax.quiver(mesh.edges["M"][:, 0],
                  mesh.edges["M"][:, 1],
                  mesh.edges["T"][:, 0],
                  mesh.edges["T"][:, 1], angles='xy', scale_units='xy', scale=5)

    if plot_normals:
        ax.quiver(mesh.edges["M"][:, 0],
                  mesh.edges["M"][:, 1],
                  mesh.edges["N"][:, 0],
                  mesh.edges["N"][:, 1], angles='xy', scale_units='xy', scale=5)

    ax.axis('equal')
    plt.show()


def _gmsh_sample_mesh():
    import pygmsh
    points = np.array([[0.0, 0.0],
                       [1.0, 0.0],
                       [1.0, 1.0],
                       [0.0, 1.0]])    

    with pygmsh.geo.Geometry() as geom:
        geom.add_polygon(points, mesh_size=0.05)
        M = geom.generate_mesh()

    return Mesh_from_meshio(M)


def _triangulation_sample_mesh():
    x = np.array([0, 1, 1, 0])
    y = np.array([0, 0, 1, 1])
    triangles = np.array([[0, 1, 2],
                          [2, 3, 0]])
    Tri = Triangulation(x=x, y=y, triangles=triangles)

    return Mesh_from_Matplotlib(Tri)


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
# def Unbounded(r = 0.1, R=5, lc=0.2) -> Mesh:
#     from pygmsh.occ import Geometry
#     with Geometry() as geom:
#         # hole = geom.add_circle( [0., 0.], r, mesh_size=lc)
#         # domain = geom.add_circle( [0., 0.], R, mesh_size=lc)
#         hole = geom.add_disk( [0., 0.], r)
#         domain = geom.add_disk( [0., 0.], R)

#         domain = geom.boolean_difference(domain, hole)
#         geom.add_physical("domain")
#         M = geom.generate_mesh()
#     return Mesh_from_meshio(M)

def _visually_test_edges(M: Mesh):
    from matplotlib.collections import LineCollection
    from matplotlib.patches import Polygon
    fig, ax = plt.subplots()
    lw = 1
    xmin, ymin = M._points.min(axis=0)
    xmax, ymax = M._points.max(axis=0)

    N_E = M.n_edges

    ax.add_collection(LineCollection(np.stack([M.edges["P"], M.edges["Q"]], axis=1), 
                                      colors='k', linewidths=lw))
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axis('equal')

    e = 0
    edge = M.edges[e]
    px, py = edge["P"]
    qx, qy = edge["Q"]
    ax.plot([px, qx], [py, qy], "b", linewidth=2*lw )
    if edge["boundary"]:
        triangle = M._triangles[edge["triangles"][0]]
        A, B, C = M._points[triangle[0]], M._points[triangle[1]], M._points[triangle[2]]
        ax.add_patch(Polygon(np.vstack([A,B,C]), facecolor='r'))
    else:
        triangle = M._triangles[edge["triangles"][0]]
        A, B, C = M._points[triangle[0]], M._points[triangle[1]], M._points[triangle[2]]
        ax.add_patch(Polygon(np.vstack([A,B,C]), facecolor='r'))

        triangle = M._triangles[edge["triangles"][1]]
        A, B, C = M._points[triangle[0]], M._points[triangle[1]], M._points[triangle[2]]
        ax.add_patch(Polygon(np.vstack([A,B,C]), facecolor='g'))


    def update_plot(e):
        ax.lines[-1].remove()
        fig.canvas.draw_idle()

        if ax.patches:
            ax.patches[-1].remove()
        if ax.patches:
            ax.patches[-1].remove()

        edge = M.edges[e]
        px, py = edge["P"]
        qx, qy = edge["Q"]
        ax.plot([px, qx], [py, qy], "b", linewidth=2*lw )
        ax.set_title(f'Edge number: {e}, boundary: {edge["boundary"]}')
        if edge["boundary"]:
            triangle = M._triangles[edge["triangles"][0]]
            A, B, C = M._points[triangle[0]], M._points[triangle[1]], M._points[triangle[2]]
            ax.add_patch(Polygon(np.vstack([A,B,C]), facecolor='r'))
        else:
            triangle = M._triangles[edge["triangles"][0]]
            A, B, C = M._points[triangle[0]], M._points[triangle[1]], M._points[triangle[2]]
            ax.add_patch(Polygon(np.vstack([A,B,C]), facecolor='r'))

            triangle = M._triangles[edge["triangles"][1]]
            A, B, C = M._points[triangle[0]], M._points[triangle[1]], M._points[triangle[2]]
            ax.add_patch(Polygon(np.vstack([A,B,C]), facecolor='g'))


        fig.canvas.draw_idle()

    def on_key(event):
        nonlocal e
        if event.key == "up":
            e = min(N_E-1, e+1)
            update_plot(e)

        elif event.key == "down":
            e = max(0, e-1)
            update_plot(e)

        elif event.key == "escape":
            plt.close(fig)

    fig.canvas.mpl_connect("key_press_event", on_key)
    update_plot(e)


    plt.show()



def Unbounded(r: float = 0.5, R: float = 5, lc: float = 0.2) -> Mesh:
    from pygmsh.geo import Geometry
    with Geometry() as geom:
        center = (0, 0)

        outer = geom.add_circle(center, R, make_surface=False, mesh_size=lc)
        inner = geom.add_circle(center, r, make_surface=False, mesh_size=lc)

        domain = geom.add_plane_surface(
            outer.curve_loop,
            holes=[inner.curve_loop]
        )

        geom.add_physical(domain, 'domain')
        geom.add_physical(outer.curve_loop, 'S_R')
        geom.add_physical(inner.curve_loop, 'S_r')

        mesh = geom.generate_mesh()

    return Mesh_from_meshio(mesh)


if __name__ == "__main__":

    # mesh = _triangulation_sample_mesh()
    # mesh = _gmsh_sample_mesh()
    # p = np.array([[0.2, 0.5],
    #               [0.5, 0.2]])
    # indexes = np.array([1, 0])

    mesh = CleanWaveGuide(lc=0.2)
    # mesh = Unbounded()

    # _triangle_id_tester(mesh)
    plot_waveguide(mesh, plot_tangents=True)
    _visually_test_edges(mesh)
