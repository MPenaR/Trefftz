import numpy as np
import pygmsh
from mesh.from_pygmsh import Mesh_from_meshio


def test_from_meshio():
    points = np.array([[0.0, 0.0],
                       [1.0, 0.0],
                       [1.0, 1.0],
                       [0.0, 1.0]])
    triangles = np.array([[0, 1, 2],
                          [2, 3, 0]])

    with pygmsh.geo.Geometry() as geom:
        for T in triangles:
            geom.add_polygon(points[T], mesh_size=10.0)
        M = geom.generate_mesh()
    mesh = Mesh_from_meshio(M)
    p = np.array([[0.2, 0.5],
                  [0.5, 0.2]])
    indexes = np.array([1, 0])
    # assert mesh.get_cell(p[0]) == indexes[0] and mesh.get_cell(p[1]) == indexes[1]
    print(mesh.get_cell(p))
    assert np.array_equal(mesh.get_cell(p), indexes)
