from trefftz.mesh import TrefftzMesh
from enum import StrEnum, IntEnum
import numpy as np

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


class Regions(IntEnum):
    OUT = 0
    BACKGROUND = 1
    OMEGA = 2

class Materials(StrEnum):
    AIR = "air"
    DIELECTRIC = "dielectric" 


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


def AnularDomain(R: float = 5., r: float = 1., lc: float = 0.3, Lc: float = 0.5, verbosity: int = 0) -> TrefftzMesh[Boundaries, Regions]:

    '''Creates an anular domain without scatterers
    It assumes the default tags for the subregions, i.e.:
    - Omega = 0
    - Sigma = 1
    '''

    geo = SplineGeometry()
    geo.AddCircle((0.0, 0.0), r, bc=Boundaries.D_OMEGA, leftdomain=Regions.OUT, rightdomain=Regions.BACKGROUND, maxh=lc)
    geo.AddCircle((0.0, 0.0), R, bc=Boundaries.SIGMA, leftdomain=Regions.BACKGROUND, rightdomain=Regions.OUT, maxh=Lc)
    ngmesh = Mesh(geo.GenerateMesh(maxh=lc, perfstepsend=verbosity))
    
    mesh = TrefftzMesh.from_ngsolve(ngmesh, boundaries=Boundaries, regions=Regions)

    return mesh


class Material():
    pass
class Dielectric(Material):
    def __init__(self, relative_permittivity: float):
        self.eps_r = relative_permittivity

class Metallic(Material):
    pass


def U(R: float = 5., lc: float = 0.3, Lc: float = 0.5, scatterer_material: Material = Metallic(), scale: float = 1., angle: float = 0., verbosity: int = 0) -> TrefftzMesh[Boundaries, Regions]:

    '''Creates an U-shaped scatterer like the one from Fresnel database:
    - Omega = 0
    - Sigma = 1
    '''

    geo = SplineGeometry()


    match scatterer_material:
        case Dielectric():
            scatterer_label = Regions.OMEGA
        case Metallic():
            scatterer_label = Regions.OUT
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
                                                                        rightdomain=Regions.BACKGROUND,
                                                                        maxh=e, bc=Boundaries.D_OMEGA)


    geo.AddCircle((0.0, 0.0), R, bc=Boundaries.SIGMA, leftdomain=Regions.BACKGROUND, rightdomain=Regions.OUT, maxh=Lc)


    ntmesh = geo.GenerateMesh(maxh=lc, perfstepsend=verbosity)

    ntmesh.SetMaterial(Regions.BACKGROUND, Materials.AIR)
    
    if isinstance(scatterer_material, Dielectric):
        ntmesh.SetMaterial(Regions.OMEGA, Materials.DIELECTRIC)

    
    mesh = TrefftzMesh.from_ngsolve(Mesh(ntmesh), boundaries=Boundaries, regions=Materials)

    return mesh