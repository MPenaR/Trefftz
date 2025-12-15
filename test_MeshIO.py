from MeshIO import Mesh
from matplotlib.tri import Triangulation
import numpy as np


def test_from_matplotlib():
    x = np.array([0, 1, 1, 0])
    y = np.array([0, 0, 1, 1])
    triangles = np.array([[0, 1, 2],
                          [2, 3, 0]])
    Tri = Triangulation(x=x, y=y, triangles=triangles)

    mesh = Mesh()
    mesh.from_matplotlib(Tri)
    assert mesh.get_cell(0.2, 0.5) == 1 and mesh.get_cell(0.5, 0.2) == 0 and mesh.n_points == 4 and mesh.n_edges == 5

