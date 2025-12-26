'''Module for importing and managing meshes'''

from typing import Final
import numpy as np
from matplotlib.tri import Triangulation
import matplotlib.pyplot as plt 
from numpy.typing import NDArray
from scipy.spatial import cKDTree
from typing import Protocol

float_array = NDArray[np.float64]
int_array = NDArray[np.int64]

# import ngsolve as ns

DIM: Final = 2

edge_dtype = [("P", np.float64, DIM), ("Q", np.float64, DIM)]


def in_triangle(P, A, B, C) -> bool:
    '''Computes if a point is inside a triangle'''
    AC = C - A
    AB = B - A
    AP = P - A

    u, v = np.linalg.solve(np.column_stack([AC,AB]), AP) #computing baricentric coordinates
    tol = 1E-16

    return (u >= -tol) and (v >= -tol) and (u + v <= 1 + tol)


class _CellLocator(Protocol):
    '''Protocol for cell locators'''
    def find_cell(self, p: float_array):
        ...


class _MatplotlibLocator(_CellLocator):
    def __init__(self, Tri: Triangulation):
        self.trifinder = Tri.get_trifinder()
    
    def find_cell(self, p: float_array):
        p_x, p_y = np.transpose(p)
        return self.trifinder(p_x, p_y)
    
class _KDTreeLocator(_CellLocator):
    def __init__(self, points, triangles):
        self.points = points
        self.triangles = triangles
        self.build_index()

    def build_index(self):
        centroids = self.points[self.triangles].mean(axis=1)
        self.tree = cKDTree(centroids)
        self.radius = np.max(np.linalg.norm(self.points[self.triangles] - centroids[:, np.newaxis,:], axis=-1))  # precomputed

    def find_cell(self, p: float_array):
        p = np.asarray(p)
        candidates = self.tree.query_ball_point(p, self.radius)


        if p.shape == (2,):
            for i in candidates:
                if in_triangle(p, *self.points[self.triangles[i]]):
                    return i
            return -1    

        elif p.ndim == 2 and p.shape[1] == 2:
            indexes = np.full(p.shape[0], dtype=np.int64, fill_value=-1)
            for j, (p_, candidates_) in enumerate(zip(p,candidates)):
                for i in candidates_:
                    if in_triangle(p_, *self.points[self.triangles[i]]):
                        indexes[j]=i
            return indexes
        else:
            raise ValueError("Input must have shape (2,) or (M, 2)")


def Mesh_from_Matplotlib(Tri: Triangulation):
        '''Returns a Mesh from a Matptlotlib.tri.Triangulation'''
        points = np.column_stack([Tri.x, Tri.y])
        edges = Tri.edges 
        triangles = Tri.triangles
        locator = _MatplotlibLocator(Tri=Tri)
        return Mesh(points=points, edges=edges, triangles=triangles, locator=locator)

def Mesh_from_meshio(mesh):
        '''Returns a Mesh from a meshio object'''
        points = mesh.points[:,0:2]
        edges = mesh.cells_dict["line"]
        triangles = mesh.cells_dict["triangle"]
        locator = _KDTreeLocator(points=points, triangles=triangles)
        cell_sets = mesh.cell_sets_dict
        return Mesh(points=points, edges=edges, triangles=triangles, locator=locator, cell_sets=cell_sets)



class Mesh():
    '''Holds only the relevant data
    as numpy structured-arrays for easy manipulation'''

    def __init__(self, points: float_array, edges: int_array, triangles: int_array, locator: _CellLocator, cell_sets : dict):
        self._points = points
        self._edges = edges
        self._triangles = triangles
        self.locator = locator
        self.construct_edge_connectivity()
        self._cell_sets = cell_sets
        self.construct_numpy_arrays()


    def construct_numpy_arrays(self):
        edges = np.zeros(self.n_edges, dtype=edge_dtype)
        points = self._points
        for e, edge in enumerate(self._edges):
            edges[e]["P"] = points[edge[0],:]
            edges[e]["Q"] = points[edge[1],:]

        self.edges = edges

        
    def construct_edge_connectivity(self):
        edges = np.array([frozenset(E) for E in self._edges])
        triangles = np.array([frozenset(T) for T in self._triangles])
        numpy_subset = np.frompyfunc(lambda A, B: A <= B, 2, 1)
        mask = numpy_subset(edges[:,None], triangles[None,:])
        self.boundary_edges_list = np.where(mask.sum(axis=1)==1)[0]
        self.inner_edges_list = np.where(mask.sum(axis=1)==2)[0]
    

    def get_cell(self, p: float_array) -> int_array:
        return self.locator.find_cell(p)
        


    @property
    def n_points(self) -> np.int64:
        return self._points.shape[0]

    @property
    def n_edges(self) -> np.int64:
        return self._edges.shape[0]
        
    def plot_mesh(self):
        plt.triplot(Triangulation(x=self._points[:,0], y=self._points[:,1], triangles=self._triangles))
        plt.show()


    def generate_Edges(self):
        pass


def _id_tester(M: Mesh):
    fig, ax = plt.subplots()

    xmin, xmax = np.min(M._points[:,0]), np.max(M._points[:,0]) 
    ymin, ymax = np.min(M._points[:,1]), np.max(M._points[:,1])

    N = 200

    x = np.linspace(xmin,xmax,N)
    y = np.linspace(ymin,ymax,N)
    X, Y = np.meshgrid(x,y)
    xy = np.column_stack([X.flatten(), Y.flatten()])
    Z = mesh.get_cell(xy).reshape(X.shape).astype(float)
    Z[Z==-1] = np.nan
    pc = ax.pcolormesh(X,Y,Z)
    lw = 2
    ax.triplot(Triangulation(x=M._points[:,0], y=M._points[:,1], triangles=M._triangles),linewidth=lw, color='k')
    for e_ID in M.boundary_edges_list:
        P, Q = M._edges[e_ID]
        p_x, p_y = M._points[P, :]
        q_x, q_y = M._points[Q, :]
        ax.plot([p_x,q_x],[p_y,q_y],'r',linewidth=lw)

    fig.colorbar(mappable=pc)
    ax.axis('equal')


    def hover_text(x, y):
        return f'ID = {M.get_cell([x, y])}'

    annot = ax.annotate(
        "",
        xy=(0, 0),
        xytext=(10, 10),
        textcoords="offset points",
        bbox=dict(boxstyle="round", fc="w"),
        arrowprops=dict(arrowstyle="->"),
    )
    annot.set_visible(False)
    annot.arrowprops = None


    def on_move(event):
        if event.inaxes != ax:
            annot.set_visible(False)
            fig.canvas.draw_idle()
            return

        xdata, ydata = event.xdata, event.ydata

        # Update annotation
        annot.xy = (xdata, ydata)
        annot.set_text(hover_text(xdata, ydata))
        annot.set_visible(True)
        fig.canvas.draw_idle()

    # Connect event
    fig.canvas.mpl_connect("motion_notify_event", on_move)

    plt.show()



def plot_waveguide(M: Mesh):
    from matplotlib.collections import LineCollection
    fig, ax = plt.subplots()

    xmin, xmax = np.min(M._points[:,0]), np.max(M._points[:,0]) 
    ymin, ymax = np.min(M._points[:,1]), np.max(M._points[:,1])

    lw = 2
    ax.triplot(Triangulation(x=M._points[:,0], y=M._points[:,1], triangles=M._triangles),linewidth=lw, color='k')

    S = M._cell_sets["S"]["line"]
    G = M._cell_sets["Gamma"]["line"]

    lines = np.stack([mesh.edges[S]["P"], mesh.edges[S]["Q"]], axis=1)
    ax.add_collection(LineCollection(lines, colors='r', linewidths=lw))

    
    # lines = [ np.vstack([mesh.edges[e_ID]["P"], mesh.edges[e_ID]["Q"]]) for e_ID in G]

    lines = np.stack([mesh.edges[G]["P"], mesh.edges[G]["Q"]], axis=1)
    ax.add_collection(LineCollection(lines, colors='b', linewidths=lw))


    ax.axis('equal')

    plt.show()




def _gmsh_sample_mesh():
    import pygmsh
    points = np.array([[0.0, 0.0],
                       [1.0, 0.0],
                       [1.0, 1.0],
                       [0.0, 1.0]])
    triangles = np.array([[0, 1, 2],
                          [2, 3, 0]])
    

    with pygmsh.geo.Geometry() as geom:
        # for T in triangles:
        #     geom.add_polygon(points[T], mesh_size=0.2)
        # M = geom.generate_mesh()
    
        square = geom.add_polygon(points, mesh_size=0.05)
        M = geom.generate_mesh()
        
    return Mesh_from_meshio(M)


def _triangulation_sample_mesh():
    x = np.array([0, 1, 1, 0])
    y = np.array([0, 0, 1, 1])
    triangles = np.array([[0, 1, 2],
                          [2, 3, 0]])
    Tri = Triangulation(x=x, y=y, triangles=triangles)

    return Mesh_from_Matplotlib(Tri)


def CleanWaveGuide(R=5, H=1, lc=0.3) -> Mesh:
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
       #  print(M.cell_sets_dict)
    return Mesh_from_meshio(M)

# def CleanWaveGuide(R=5, H=1, lc=0.2) -> Mesh:
#     from pygmsh.geo import Geometry
#     points = np.array([[-R, 0.],
#                        [ R, 0.],
#                        [ R, H ],
#                        [-R, H ]])
#     with Geometry() as geom:
#         geom.add_polygon(points, mesh_size=lc)
#         M = geom.generate_mesh()
#     return Mesh_from_meshio(M)




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

def Unbounded(r = 0.5, R=5, lc=0.2) -> Mesh:
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
    mesh = _gmsh_sample_mesh()
    p = np.array([[0.2, 0.5],
                  [0.5, 0.2]])
    indexes = np.array([1, 0])

    mesh = CleanWaveGuide()
    # mesh = Unbounded()
    

    # _id_tester(mesh)
    plot_waveguide(mesh)