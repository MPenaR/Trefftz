from trefftz.mesh import TrefftzMesh
from enum import StrEnum

try:
    from netgen.geom2d import SplineGeometry
    from ngsolve import Mesh
    
except ImportError as e:
    raise RuntimeError(
        "NGSolve is required for this feature. "
        "Install it with `pip install trefftz[ngsolve]` or `uv pip install trefftz[ngsolve]`."
    ) from e


# class WaveguideRegions(StrEnum):
#     OMEGA = "Omega"
#     GAMMA = "Gamma"
#     SIGMA_L = "Sigma_L"
#     SIGMA_R = "Sigma_R"

class BoundaryRegions(StrEnum):
    GAMMA = "Gamma"
    SIGMA_L = "Sigma_L"
    SIGMA_R = "Sigma_R"


def Empty(H: float = 1., R: float = 5., lc: float = 0.3, verbosity: int = 0) -> TrefftzMesh[BoundaryRegions]:

    '''Creates a domain corresponging to a waveguide without scatterers.
    It assumes the default tags for the subregions, i.e.:
    - Gamma = "Gamma"
    - Sigma_L = "Sigma_L"
    - Sigma_R = "Sigma_R"
    '''

    geo = SplineGeometry()
    p0 = geo.AddPoint(-R, 0.)
    p1 = geo.AddPoint( R, 0.)
    p2 = geo.AddPoint( R, H)
    p3 = geo.AddPoint(-R, H)
 
    bottom = geo.Append(["line", p0, p1], bc=BoundaryRegions.GAMMA)
    right  = geo.Append(["line", p1, p2], bc=BoundaryRegions.SIGMA_R)
    top    = geo.Append(["line", p2, p3], bc=BoundaryRegions.GAMMA)
    left   = geo.Append(["line", p3, p0], bc=BoundaryRegions.SIGMA_L)

    ngmesh = Mesh(geo.GenerateMesh(maxh=lc, perfstepsend=verbosity))
    
    mesh = TrefftzMesh.from_ngsolve(ngmesh, boundary_regions=BoundaryRegions)

    return mesh



