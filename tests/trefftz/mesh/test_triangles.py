from trefftz.mesh.from_pygmsh import Mesh_from_meshio
import numpy as np

def test_triangle_area():
    from pygmsh.geo import Geometry
    with Geometry() as geom:
        geom.add_circle((0, 0), 1., mesh_size=0.05)
        mesh = Mesh_from_meshio(geom.generate_mesh())
    assert np.isclose(np.sum(mesh.triangles["area"]), np.pi, rtol=1e-2)