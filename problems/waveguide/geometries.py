from typing import TYPE_CHECKING


from trefftz.mesh import TrefftzMesh
from enum import IntEnum


from trefftz.mesh.from_pygmsh import Mesh_from_meshio
from pygmsh.geo import Geometry

import gmsh


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


def CleanWaveGuide2(R: float = 5, H: float = 1, lc: float = 0.3) -> "TrefftzMesh":
    '''Constructs a waveguide mesh without scatterers'''
    gmsh.initialize()
    gmsh.model.add("Waveguide")
    p0 = gmsh.model.geo.addPoint(-R, 0., 0., lc)
    p1 = gmsh.model.geo.addPoint( R, 0., 0., lc)
    p2 = gmsh.model.geo.addPoint( R,  H, 0., lc)
    p3 = gmsh.model.geo.addPoint(-R,  H, 0., lc)
    
    bottom = gmsh.model.geo.addLine(p0, p1)
    right  = gmsh.model.geo.addLine(p1, p2)
    top    = gmsh.model.geo.addLine(p2, p3)
    left   = gmsh.model.geo.addLine(p3, p0)

    boundary = gmsh.model.geo.addCurveLoop([bottom, right, top, left])
    domain = gmsh.model.geo.addPlaneSurface([boundary])
    gmsh.model.addPhysicalGroup(2, [domain], 0, "Omega")
    gmsh.model.addPhysicalGroup(1, [bottom, top], 1, "Gamma")
    gmsh.model.addPhysicalGroup(1, [left], 3, "Sigma_L")
    gmsh.model.addPhysicalGroup(1, [right], 4, "Sigma_R")
    gmsh.model.addPhysicalGroup(1, [left, right], 2, "Sigma")

    
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    gmsh.fltk.run()
    gmsh.write('CleanWaveguide.msh')
    gmsh.finalize()
    return TrefftzMesh.from_gmsh('CleanWaveguide.msh')
