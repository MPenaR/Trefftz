from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from trefftz.mesh import TrefftzMesh


from trefftz.mesh.from_pygmsh import Mesh_from_meshio
from pygmsh.geo import Geometry


def CleanWaveGuide(R: float = 5, H: float = 1, lc: float = 0.3) -> "TrefftzMesh":
    '''Constructs a waveguide mesh without scatterers'''
    with Geometry() as geom:
        p0 = geom.add_point([-R, 0.], mesh_size=lc)
        p1 = geom.add_point([ R, 0.], mesh_size=lc)
        p2 = geom.add_point([ R, H ], mesh_size=lc)
        p3 = geom.add_point([-R, H ], mesh_size=lc)

        bottom = geom.add_line( p0, p1)
        right  = geom.add_line(p1, p2)
        top    = geom.add_line(p2, p3)
        left   = geom.add_line(p3, p0)

        boundary = geom.add_curve_loop([bottom, right, top, left])
        domain = geom.add_plane_surface(boundary)

        geom.add_physical(domain, "Omega")
        geom.add_physical([bottom, top], "Gamma")
        geom.add_physical(left, "S_L")
        geom.add_physical(right, "S_R")
        geom.add_physical([left, right], "S")

        M = geom.generate_mesh()
    return Mesh_from_meshio(M)
