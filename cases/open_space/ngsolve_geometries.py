from trefftz.mesh import TrefftzMesh
from enum import StrEnum, IntEnum
import numpy as np
from trefftz.dg.materials import Dielectric, Metallic, Material

try:
    from netgen.geom2d import SplineGeometry
    from ngsolve import Mesh
    
except ImportError as e:
    raise RuntimeError(
        "NGSolve is required for this feature. "
        "Install it with `pip install trefftz[ngsolve]` or `uv pip install trefftz[ngsolve]`."
    ) from e


class Boundaries(StrEnum):
    SIGMA = "Sigma"
    D_OMEGA = "d_Omega"


class Labels(IntEnum):
    OUT = 0
    BACKGROUND = 1
    SCATTERER = 2

class Regions(StrEnum):
    BACKGROUND = "Background"
    OMEGA = "Omega" 


def CleanCircle(R: float = 5., lc: float = 0.3, verbosity: int = 0) -> TrefftzMesh[Boundaries, Regions]:

    '''Creates a circular domain without scatterers
    It assumes the default tags for the subregions, i.e.:
    - Omega = 0
    - Sigma = 1
    '''

    geo = SplineGeometry()
    geo.AddCircle((0.0, 0.0), R, bc=Boundaries.SIGMA, maxh=lc)

    ngmesh = Mesh(geo.GenerateMesh(maxh=lc, perfstepsend=verbosity))
    
    mesh = TrefftzMesh.from_ngsolve(ngmesh, boundaries=Boundaries, regions=Regions)

    return mesh


def AnularDomain(R: float = 5., r: float = 1., scatterer_material: Material = Metallic(), lc: float = 0.3, Lc: float = 0.5, verbosity: int = 0) -> TrefftzMesh[Boundaries, Regions]:

    '''Creates an anular domain without scatterers
    It assumes the default tags for the subregions, i.e.:
    - Omega = 0
    - Sigma = 1
    '''

    match scatterer_material:
        case Dielectric():
            scatterer_label = Labels.SCATTERER
        case Metallic():
            scatterer_label = Labels.OUT
        case _:
            raise TypeError(f"Unsupported scatterer material: {type(scatterer_material).__name__}")

    geo = SplineGeometry()
    geo.AddCircle((0.0, 0.0), r, bc=Boundaries.D_OMEGA, leftdomain=scatterer_label, rightdomain=Labels.BACKGROUND, maxh=lc)
    geo.AddCircle((0.0, 0.0), R, bc=Boundaries.SIGMA, leftdomain=Labels.BACKGROUND, rightdomain=Labels.OUT, maxh=Lc)
    
    ntmesh = geo.GenerateMesh(maxh=lc, perfstepsend=verbosity)

    ntmesh.SetMaterial(Labels.BACKGROUND, Regions.BACKGROUND)
    
    if isinstance(scatterer_material, Dielectric):
        ntmesh.SetMaterial(Labels.SCATTERER, Regions.OMEGA)

    mesh = TrefftzMesh.from_ngsolve(Mesh(ntmesh), boundaries=Boundaries, regions=Regions)

    return mesh

def U(R: float = 5., lc: float = 0.3, Lc: float = 0.5, scatterer_material: Material = Metallic(), scale: float = 1., angle: float = 0., verbosity: int = 0) -> TrefftzMesh[Boundaries, Regions]:

    '''Creates an U-shaped scatterer like the one from Fresnel database:
    - Omega = 0
    - Sigma = 1
    '''

    geo = SplineGeometry()


    match scatterer_material:
        case Dielectric():
            scatterer_label = Labels.SCATTERER
        case Metallic():
            scatterer_label = Labels.OUT
        case _:
            raise TypeError(f"Unsupported scatterer material: {type(scatterer_material).__name__}")



    # scatterer
    L = 0.080
    H = 0.050
    e = 0.005
    vertices = np.array([[ -L/2    , -H/2    ],
                         [  L/2    , -H/2    ],
                         [  L/2    , -H/2 + e],
                         [ -L/2 + e, -H/2 + e],
                         [ -L/2 + e,  H/2 - e],
                         [  L/2    ,  H/2 - e],
                         [  L/2    ,  H/2    ],
                         [ -L/2    ,  H/2    ]])

    Q = np.array( [[np.cos(angle), -np.sin(angle)],
                   [np.sin(angle),  np.cos(angle)]])
    
    vertices = vertices*scale  # scaling the U
    vertices = np.matmul(vertices, Q.T)   # rotating the U

    corners = [ geo.AppendPoint(*v) for v in vertices]
    

    for i in range(len(corners)):
        geo.Append(["line", corners[i], corners[(i+1) % len(corners)]], leftdomain=scatterer_label,
                                                                        rightdomain=Labels.BACKGROUND,
                                                                        maxh=e, bc=Boundaries.D_OMEGA)


    geo.AddCircle((0.0, 0.0), R, bc=Boundaries.SIGMA, leftdomain=Labels.BACKGROUND, rightdomain=Labels.OUT, maxh=Lc)


    ntmesh = geo.GenerateMesh(maxh=lc, perfstepsend=verbosity)

    ntmesh.SetMaterial(Labels.BACKGROUND, Regions.BACKGROUND)
    
    if isinstance(scatterer_material, Dielectric):
        ntmesh.SetMaterial(Labels.SCATTERER, Regions.OMEGA)

    
    mesh = TrefftzMesh.from_ngsolve(Mesh(ntmesh), boundaries=Boundaries, regions=Regions)

    return mesh
        