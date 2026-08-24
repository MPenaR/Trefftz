from trefftz.mesh import TrefftzMesh
from enum import StrEnum
import numpy as np

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





def U(R: float = 5., lc: float = 0.3, Lc: float = 0.5, scale: float = 1., angle: float = 0., verbosity: int = 0) -> TrefftzMesh[Regions]:

    '''Creates an U-shaped scatterer like the one from Fresnel database:
    - Omega = 0
    - Sigma = 1
    '''

    geo = SplineGeometry()


    # match scatterer_material:
    #     case Dielectric():
    #         scatterer_label = Region.SCATTERER
    #     case Metallic():
    #         scatterer_label = Region.OUT
    #     case _:
    #         raise TypeError(f"Unsupported scatterer material: {type(scatterer_material).__name__}")


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
    
    # for i in range(len(corners)):
    #     geo.Append(["line", corners[i], corners[(i+1) % len(corners)]], leftdomain=scatterer_label, rightdomain=Region.BACKGROUND, maxh=e, bc=Regions.d_Omega)

    for i in range(len(corners)):
        geo.Append(["line", corners[i], corners[(i+1) % len(corners)]], leftdomain=0, rightdomain=1, maxh=e, bc=Regions.d_Omega)



    # mesh.SetMaterial(Region.BACKGROUND, materials.AIR)
    # if isinstance(scatterer_material, Dielectric):
    #     mesh.SetMaterial(Region.SCATTERER, materials.DIELECTRIC)

    geo.AddCircle((0.0, 0.0), R, bc=Regions.SIGMA, leftdomain=1, rightdomain=0, maxh=Lc)


    ngmesh = Mesh(geo.GenerateMesh(maxh=lc, perfstepsend=verbosity))
    
    mesh = TrefftzMesh.from_ngsolve(ngmesh, boundary_regions=Regions)

    return mesh