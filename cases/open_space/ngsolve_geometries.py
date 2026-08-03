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


class Regions(StrEnum):
    OMEGA = "Omega"
    SIGMA = "Sigma"
    d_Omega = "d_Omega"


def CleanCircle(R: float = 5., lc: float = 0.3, verbosity: int = 0) -> TrefftzMesh[Regions]:

    '''Creates a circular domain without scatterers
    It assumes the default tags for the subregions, i.e.:
    - Omega = 0
    - Sigma = 1
    '''

    geo = SplineGeometry()
    geo.AddCircle((0.0, 0.0), R, bc=Regions.SIGMA, maxh=lc)

    ngmesh = Mesh(geo.GenerateMesh(maxh=lc, perfstepsend=verbosity))
    
    mesh = TrefftzMesh.from_ngsolve(ngmesh, boundary_regions=Regions)

    return mesh


def AnularDomain(R: float = 5., r: float = 1., lc: float = 0.3, Lc: float = 0.5, verbosity: int = 0) -> TrefftzMesh[Regions]:

    '''Creates an anular domain without scatterers
    It assumes the default tags for the subregions, i.e.:
    - Omega = 0
    - Sigma = 1
    '''

    geo = SplineGeometry()
    geo.AddCircle((0.0, 0.0), r, bc=Regions.d_Omega, leftdomain=0, rightdomain=1, maxh=lc)
    geo.AddCircle((0.0, 0.0), R, bc=Regions.SIGMA, leftdomain=1, rightdomain=0, maxh=Lc)
    ngmesh = Mesh(geo.GenerateMesh(maxh=lc, perfstepsend=verbosity))
    
    mesh = TrefftzMesh.from_ngsolve(ngmesh, boundary_regions=Regions)

    return mesh