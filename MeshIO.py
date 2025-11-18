'''Module for importing and managing meshes'''

from typing import Final
import numpy as np
from matplotlib.tri import Triangulation
import matplotlib.pyplot as plt 

DIM: Final = 2


class Mesh():
    '''Holds only the relevant data
    as numpy structured-arrays for easy manipulation'''

    def from_matplotlib(self, Tri: Triangulation ):
        self.points = np.column_stack([Tri.x, Tri.y])
        self.get_cell = Tri.get_trifinder()





if __name__ == "__main__":
    x = np.array([0, 1, 1, 0])
    y = np.array([0, 0, 1, 1])
    triangles = np.array([[0, 1, 2],
                          [2, 3, 0]])
    Tri = Triangulation(x=x, y=y, triangles=triangles)

    mesh = Mesh()
    mesh.from_matplotlib(Tri)
    print(Tri.edges)

    print(mesh.points)

    plt.triplot(Tri)
    plt.show()   

