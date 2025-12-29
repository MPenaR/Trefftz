import numpy as np
import pygmsh
from trefftz.mesh.from_pygmsh import Mesh_from_meshio


# def gmsh_sample_mesh():
#     from pygmsh.geo import Geometry
#     points = np.array([[0.0, 0.0],
#                        [1.0, 0.0],
#                        [1.0, 1.0],
#                        [0.0, 1.0]])    

#     with Geometry() as geom:
#         geom.add_polygon(points, mesh_size=0.05)
#         M = geom.generate_mesh()

#     return Mesh_from_meshio(M)


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
    assert np.array_equal(mesh.get_cell(p), indexes)
