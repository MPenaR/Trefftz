'''Module for importing and managing meshes'''

from typing import Final
import numpy as np
from numpy.linalg import norm
from matplotlib.tri import Triangulation
import matplotlib.pyplot as plt 
from typing import Protocol

from numpy.typing import NDArray
float_array = NDArray[np.float64]
int_array = NDArray[np.int64]


DIM: Final = 2

edge_dtype = [("P", np.float64, DIM),
              ("Q", np.float64, DIM),
              ("T", np.float64, DIM),
              ("N", np.float64, DIM),
              ("M", np.float64, DIM),
              ("l", float),
              ("boundary", bool),
              ("triangles", np.int32, 2)]


class CellLocator(Protocol):
    '''Protocol for cell locators'''
    def find_cell(self, p: float_array) -> int_array | int:
        ...


class Mesh():
    '''Holds only the relevant data
    as numpy structured-arrays for easy manipulation'''

    def __init__(self, points: float_array, edges: int_array, triangles: int_array,
                 edge2triangles: int_array,
                 locator: CellLocator, cell_sets: dict[str, dict[str, int_array]]):
        self._points = points
        self._edges = edges
        self._triangles = triangles
        self.locator = locator
        self._cell_sets = cell_sets
        self._edge2triangles = edge2triangles
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
        edges["triangles"] = self._edge2triangles
        edges["boundary"] = edges["triangles"][:, 1] == -1
        self.edges = edges

    def get_cell(self, p: float_array) -> int_array | int:
        return self.locator.find_cell(p)

    @property
    def n_points(self) -> int:
        return self._points.shape[0]

    @property
    def n_edges(self) -> int:
        return self._edges.shape[0]



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


