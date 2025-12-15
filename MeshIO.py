'''Module for importing and managing meshes'''

from typing import Final
import numpy as np
from matplotlib.tri import Triangulation
import matplotlib.pyplot as plt 

# import ngsolve as ns

DIM: Final = 2


class Mesh():
    '''Holds only the relevant data
    as numpy structured-arrays for easy manipulation'''

    def from_matplotlib(self, Tri: Triangulation):
        '''Constructs the object from a Matptlotlib.tri.Triangulation'''
        self._points = np.column_stack([Tri.x, Tri.y])
        self._triangles = Tri.triangles
        self._edges = Tri.edges 
        self.get_cell = Tri.get_trifinder()
    
    # def from_netgen(self, Tri: ns.Mesh):
    #     pass

    @property
    def n_points(self) -> np.int64:
        return self._points.shape[0]

    @property
    def n_edges(self) -> np.int64:
        return self._edges.shape[0]
        

    def generate_Edges(self):
        pass





if __name__ == "__main__":
    x = np.array([0, 1, 1, 0])
    y = np.array([0, 0, 1, 1])
    triangles = np.array([[0, 1, 2],
                          [2, 3, 0]])
    Tri = Triangulation(x=x, y=y, triangles=triangles)

    mesh = Mesh()
    mesh.from_matplotlib(Tri)
    
    plt.triplot(Tri)
    plt.show()   

