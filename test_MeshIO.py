from MeshIO import Mesh, Mesh_from_Matplotlib
from matplotlib.tri import Triangulation
import numpy as np


def test_from_matplotlib():
    x = np.array([0, 1, 1, 0])
    y = np.array([0, 0, 1, 1])
    triangles = np.array([[0, 1, 2],
                          [2, 3, 0]])
    Tri = Triangulation(x=x, y=y, triangles=triangles)

    mesh = Mesh_from_Matplotlib(Tri=Tri)
    # assert mesh.get_cell(np.array([0.2, 0.5])) == 1 and mesh.get_cell(np.array([0.5, 0.2])) == 0 and mesh.n_points == 4 and mesh.n_edges == 5
    p = np.array([[0.2, 0.5],
                  [0.5, 0.2]])
    indexes = np.array([1, 0])
    assert np.array_equal(mesh.get_cell(p), indexes) and mesh.n_points == 4 and mesh.n_edges == 5
