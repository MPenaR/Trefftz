from trefftz.mesh import TrefftzMesh
from enum import IntEnum

try:
    import gmsh
except ImportError as e:
    raise RuntimeError(
        "gmsh is required for this feature. "
        "Install it with `pip install trefftz[gmsh]` or `uv pip install trefftz[gmsh]`."
    ) from e


class WaveguideRegions(IntEnum):
    OMEGA = 0
    GAMMA = 1
    SIGMA_L = 2
    SIGMA_R = 3


def CleanWaveguide(H: float = 1., R: float = 5., lc: float = 0.3, verbose: bool = False) -> TrefftzMesh[WaveguideRegions]:
    '''Creates a domain corresponging to a waveguide without scatterers.
    It assumes the default tags for the subregions, i.e.:
    - Omega = 0
    - Gamma = 1
    - Sigma = 2
    '''

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", int(verbose))
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
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(2, [domain], WaveguideRegions.OMEGA, "Omega")
    gmsh.model.addPhysicalGroup(1, [bottom, top], WaveguideRegions.GAMMA, "Gamma")
    gmsh.model.addPhysicalGroup(1, [left], WaveguideRegions.SIGMA_L, "Sigma_L")
    gmsh.model.addPhysicalGroup(1, [right], WaveguideRegions.SIGMA_R, "Sigma_R")

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    mesh = TrefftzMesh.from_gmsh(gmsh.model, boundary_regions=WaveguideRegions)  # gmsh needs to be initialized for this to work

    gmsh.finalize()

    return mesh
