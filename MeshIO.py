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
        self.radius = np.max(np.linalg.norm(self.points[self._triangles] - centroids[:, np.newaxis,:], axis=-1))  # precomputed

    def find_cell(self, p: float_array):
        candidates = self.tree.query_ball_point(p, self.radius)
        for i in candidates:
            if in_triangle(p, *self.points[self.triangles[i]]):
                return i
        return -1
    
    def find_cell(self, p: float_array):
        p_x, p_y = np.transpose(p)
        return self.trifinder(p_x, p_y)



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
        locator = _GeneralLocator(points=points, triangles=triangles)
        return Mesh(points=points, edges=edges, triangles=triangles, locator=locator)



class Mesh():
    '''Holds only the relevant data
    as numpy structured-arrays for easy manipulation'''

    def __init__(self, points: float_array, edges: int_array, triangles: int_array, locator: _CellLocator):
        self._points = points
        self._edges = edges
        self._triangles = triangles
        self.locator = locator
    

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
        

    def generate_Edges(self):
        pass





if __name__ == "__main__":
    # x = np.array([0, 1, 1, 0])
    # y = np.array([0, 0, 1, 1])
    # triangles = np.array([[0, 1, 2],
    #                       [2, 3, 0]])
    # Tri = Triangulation(x=x, y=y, triangles=triangles)

    # mesh = Mesh()
    # mesh.from_matplotlib(Tri)
    import pygmsh
    
    with pygmsh.geo.Geometry() as geom:
        geom.add_polygon(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [1.0, 1.0],
                [0.0, 1.0],
            ],
            mesh_size=1.0,
        )
        mesh = geom.generate_mesh()
    M = Mesh()
    M.from_gmsh(mesh)
    x = M._points[:,0]
    y = M._points[:,1]
    triangles = M._triangles
    Tri = Triangulation(x=x,y=y,triangles=triangles)
    plt.triplot(Tri)
    plt.show()   

    for T in M._triangles:
        print(in_triangle(np.array([0.5,0.25]), *(M._points[T])))

