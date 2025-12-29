import numpy as np
from matplotlib.tri import Triangulation
from mesh.from_matplotlib import Mesh_from_Matplotlib
import pytest


# def triangulation_sample_mesh():
#     x = np.array([0, 1, 1, 0])
#     y = np.array([0, 0, 1, 1])
#     triangles = np.array([[0, 1, 2],
#                           [2, 3, 0]])
#     Tri = Triangulation(x=x, y=y, triangles=triangles)

#     return Mesh_from_Matplotlib(Tri)



@pytest.mark.xfail
def test_from_matplotlib():
    x = np.array([0, 1, 1, 0])
    y = np.array([0, 0, 1, 1])
    triangles = np.array([[0, 1, 2],
                          [2, 3, 0]])
    Tri = Triangulation(x=x, y=y, triangles=triangles)

    mesh = Mesh_from_Matplotlib(Tri=Tri)
    p = np.array([[0.2, 0.5],
                  [0.5, 0.2]])
    indexes = np.array([1, 0])
    assert np.array_equal(mesh.get_cell(p), indexes) and mesh.n_points == 4 and mesh.n_edges == 5
