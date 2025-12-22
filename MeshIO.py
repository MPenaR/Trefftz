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
        edges = mesh.cells[1].data
        triangles = mesh.cells[1].data
        locator = _KDTreeLocator(points=points, triangles=triangles)
        return Mesh(points=points, edges=edges, triangles=triangles, locator=locator)



class Mesh():
    '''Holds only the relevant data
    as numpy structured-arrays for easy manipulation'''

    def __init__(self, points: float_array, edges: int_array, triangles: int_array, locator: _CellLocator):
        self._points = points
        self._edges = edges
        self._triangles = triangles
        self.locator = locator


    def construct_edge_connectivity(self):
        mask = np.isin(self._edges[:, None, :], self._triangles[None, :, :])  # (NE, NT, 2)
        belongs = mask.all(axis=2)  # (NE, NT)
        self.boundary_edges_list = np.where(belongs.sum(axis=1)==1)[0]
        self.inner_edges_list = np.where(belongs.sum(axis=1)==2)[0]
        

    

    def get_cell(self, p: float_array):
        return self.locator.find_cell(p)

    # def from_matplotlib(self, Tri: Triangulation):
    #     '''Constructs the object from a Matptlotlib.tri.Triangulation'''
    #     self._points = np.column_stack([Tri.x, Tri.y])
    #     self._edges = Tri.edges 
    #     self._triangles = Tri.triangles
    #     self.get_cell = Tri.get_trifinder()
    
    # # def from_netgen(self, Tri: ns.Mesh):
    # #     pass
    # def from_gmsh(self, mesh):
    #     self._points = mesh.points[:,0:2]
    #     mesh._edges = mesh.cells[1].data
    #     self._triangles = mesh.cells[1].data
    #     self._centroids = np.mean(self._points[self._triangles], axis=1)



    #     def get_cell(P: float_array) -> int_array:                
    #         tree = cKDTree(self._centroids)
    #         neighbours_radius = np.max(np.linalg.norm(self._points[self._triangles] - self._centroids[:, np.newaxis,:], axis=-1))
        


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

    x = np.linspace(0,1,100)
    y = x
    X, Y = np.meshgrid(x,y)
    xy = np.column_stack([X.flatten(), Y.flatten()])
    Z = mesh.get_cell(xy).reshape(X.shape)
    pc = ax.pcolormesh(X,Y,Z)
    ax.triplot(Triangulation(x=M._points[:,0], y=M._points[:,1], triangles=M._triangles), color='k')
    fig.colorbar(mappable=pc)


    # Your function: f(x, y) -> text
    def hover_text(x, y):
        return f'ID = {M.get_cell([x, y])}'

    # Create annotation (label)
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





if __name__ == "__main__":
    # x = np.array([0, 1, 1, 0])
    # y = np.array([0, 0, 1, 1])
    # triangles = np.array([[0, 1, 2],
    #                       [2, 3, 0]])
    # Tri = Triangulation(x=x, y=y, triangles=triangles)

    # mesh = Mesh_from_Matplotlib(Tri)
    # x = np.linspace(0,1)
    # y = x
    # X, Y = np.meshgrid(x,y)
    # xy = np.column_stack([X.flatten(), Y.flatten()])
    # Z = mesh.get_cell(xy).reshape(X.shape)
    # plt.pcolormesh(X,Y,Z)
    # plt.show()

    import pygmsh
    
    points = np.array([[0.0, 0.0],
                       [1.0, 0.0],
                       [1.0, 1.0],
                       [0.0, 1.0]])
    triangles = np.array([[0, 1, 2],
                          [2, 3, 0]])

    with pygmsh.geo.Geometry() as geom:
        for T in triangles:
            geom.add_polygon(points[T], mesh_size=0.2)
        M = geom.generate_mesh()
    mesh = Mesh_from_meshio(M)
    p = np.array([[0.2, 0.5],
                  [0.5, 0.2]])
    indexes = np.array([1, 0])
    print(mesh.get_cell(p[0]))
    print(mesh.get_cell([0.5,0.2]))
    # for T in M._triangles:
    #     print(in_triangle(np.array([0.5,0.25]), *(M._points[T])))

    _id_tester(mesh)