'''Module for defining a waveguide class, as it will be, I think, the most usefull'''
from trefftz.mesh import TrefftzMesh, EdgeType
# from trefftz.mesh.from_pygmsh import Mesh_from_meshio
from typing import Optional, Callable, Any
# from pygmsh.geo import Geometry
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from regions import Region
from problems.base import Domain
import gmsh

from trefftz.mesh.readers import GmshArrays
from problems import Problem
from trefftz.dg.fluxes import FluxType
from problems.waveguide.exact_solutions import Mode
from trefftz.numpy_types import float_array, complex_array

# from scipy.sparse import csc_matrix
# from trefftz.numpy_types import float_array




# class FluxType(IntEnum):  # should go in flux types or in fluxes
#     TRANSMISSION = 0
#     SOUNDHARD = 1
#     SOUNDSOFT = 2
#     RADIATION = 3


class WaveguideDomain(Domain):
    def __init__(self, mesh: "TrefftzMesh", regions: Region, R: float = 5., H: float = 1.):
        self.R = R
        self.H = H
        self.mesh = mesh
        self.regions = regions




class Scatterer:
    pass

# class Waveguide:
#     def __init__(self, H: float = 1., R: float = 5., lc: float = 0.3,
#                  scatterers: Optional[tuple[Scatterer]] = None):
#         self.R = R
#         self.H = H
#         with Geometry() as geom:
#             p0 = geom.add_point([-R, 0.], mesh_size=lc)
#             p1 = geom.add_point([ R, 0.], mesh_size=lc)
#             p2 = geom.add_point([ R, H ], mesh_size=lc)
#             p3 = geom.add_point([-R, H ], mesh_size=lc)

#             bottom = geom.add_line( p0, p1)
#             right  = geom.add_line(p1, p2)
#             top    = geom.add_line(p2, p3)
#             left   = geom.add_line(p3, p0)

#             boundary = geom.add_curve_loop([bottom, right, top, left])

#             if scatterers:
#                 pass


#             domain = geom.add_plane_surface(boundary)

#             geom.add_physical(domain, label="Omega")
#             geom.add_physical([bottom, top], label="Gamma")
#             geom.add_physical(left, label="S_L")
#             geom.add_physical(right, label="S_R")
#             geom.add_physical([left, right], label="S")
#             self._domain = Mesh_from_meshio(geom.generate_mesh())
    
#     @property
#     def domain(self) -> "TrefftzMesh":
#         return self._domain
    

#     def plot(self, figsize: Optional[tuple[int, int]] = (16, 2), line_width: Optional[int] = 1):
#         from matplotlib.collections import LineCollection
#         _, ax = plt.subplots(figsize=figsize)
#         # ax.triplot(Triangulation(x=M._points[:,0], y=M._points[:,1], triangles=M._triangles),linewidth=lw, color='k')

#         S = self.domain._cell_sets["S"]["line"]
#         G = self.domain._cell_sets["Gamma"]["line"]  # it allows for multidimensional subsets

#         inner = np.where(self.domain.edges["type"] == EdgeType.INNER)[0]

#         ax.add_collection(LineCollection(np.stack([self.domain.edges[inner]["P"],
#                                                    self.domain.edges[inner]["Q"]], axis=1),
#                                                    colors='k', linewidths=line_width))

#         ax.add_collection(LineCollection(np.stack([self.domain.edges[S]["P"],
#                                                    self.domain.edges[S]["Q"]], axis=1),
#                                                    colors='r', linewidths=line_width))

#         ax.add_collection(LineCollection(np.stack([self.domain.edges[G]["P"],
#                                                    self.domain.edges[G]["Q"]], axis=1),
#                                                    colors='b', linewidths=line_width))


#         ax.axis('equal')
#         ax.axis('off')
#         plt.show()
    

#     def plot_field(self, u: Callable[[NDArray[Any], NDArray[Any]], NDArray[Any]],
#                    N: int = 100, figsize: Optional[tuple[int, int]] = (16, 2), real_part: bool = False):
#         x = np.linspace(-self.R, self.R, N)
#         y = np.linspace(0., self.H, N)
#         X, Y = np.meshgrid(x, y)
#         Z = u(X, Y)

#         if real_part:
#             Z = np.real(Z)

#         _, ax = plt.subplots(figsize=figsize)

#         ax.pcolorfast((-self.R, self.R), (0., self.H), Z)
#         plt.show()



# class Waveguide:
#     def __init__(self, H: float = 1., R: float = 5., lc: float = 0.3,
#                  scatterers: Optional[tuple[Scatterer]] = None, verbose: bool = False):
#         self.R = R
#         self.H = H
#         self.verbose = verbose
#         gmsh.initialize()
#         gmsh.option.setNumber("General.Terminal", int(verbose))
#         gmsh.model.add("Waveguide")
#         p0 = gmsh.model.geo.addPoint(-R, 0., 0., lc)
#         p1 = gmsh.model.geo.addPoint( R, 0., 0., lc)
#         p2 = gmsh.model.geo.addPoint( R,  H, 0., lc)
#         p3 = gmsh.model.geo.addPoint(-R,  H, 0., lc)
        
#         bottom = gmsh.model.geo.addLine(p0, p1)
#         right  = gmsh.model.geo.addLine(p1, p2)
#         top    = gmsh.model.geo.addLine(p2, p3)
#         left   = gmsh.model.geo.addLine(p3, p0)

#         boundary = gmsh.model.geo.addCurveLoop([bottom, right, top, left])
#         domain = gmsh.model.geo.addPlaneSurface([boundary])
#         gmsh.model.addPhysicalGroup(2, [domain], Region.OMEGA, "Omega")
#         gmsh.model.addPhysicalGroup(1, [bottom, top], Region.GAMMA, "Gamma")
#         # gmsh.model.addPhysicalGroup(1, [left], Region.SIGMA_L, "Sigma_L")
#         # gmsh.model.addPhysicalGroup(1, [right], Region.SIGMA_R, "Sigma_R")
#         gmsh.model.addPhysicalGroup(1, [left, right], Region.SIGMA, "Sigma")
        
#         gmsh.model.geo.synchronize()
#         gmsh.model.mesh.generate(2)
#         # gmsh.fltk.run()
#         # gmsh.write('CleanWaveguide.msh')
#         # gmsh.finalize()
#         # self._domain = TrefftzMesh.from_gmsh('CleanWaveguide.msh')
#         points, edges, triangles, edges2triangles, locator, cell_sets = GmshArrays(gmsh.model)
#         gmsh.finalize()

#         # construction cell sets (regions)

#         self._domain = TrefftzMesh(points, edges, triangles, edges2triangles, locator, cell_sets)
#         self._regions = Region
#         return
    
#     @property
#     def domain(self) -> "TrefftzMesh":
#         return self._domain
    

#     def plot(self, figsize: Optional[tuple[int, int]] = (16, 2), line_width: Optional[int] = 1):
#         from matplotlib.collections import LineCollection
#         _, ax = plt.subplots(figsize=figsize)
#         # ax.triplot(Triangulation(x=M._points[:,0], y=M._points[:,1], triangles=M._triangles),linewidth=lw, color='k')
    
#         # S = self.domain._cell_sets["S"]["line"]
#         # G = self.domain._cell_sets["Gamma"]["line"]  # it allows for multidimensional subsets
#         # inner = np.where(self.domain.edges["type"] == EdgeType.INNER)[0]

#         S = np.where(self.domain.edges["region"] == Region.SIGMA)[0]
#         G = np.where(self.domain.edges["region"] == Region.GAMMA)[0]  
#         inner = np.where(self.domain.edges["type"] == EdgeType.INNER)[0]
    


#         ax.add_collection(LineCollection(np.stack([self.domain.edges[inner]["P"],
#                                                    self.domain.edges[inner]["Q"]], axis=1),
#                                                    colors='k', linewidths=line_width))

#         ax.add_collection(LineCollection(np.stack([self.domain.edges[S]["P"],
#                                                    self.domain.edges[S]["Q"]], axis=1),
#                                                    colors='r', linewidths=line_width))

#         ax.add_collection(LineCollection(np.stack([self.domain.edges[G]["P"],
#                                                    self.domain.edges[G]["Q"]], axis=1),
#                                                    colors='b', linewidths=line_width))


#         ax.axis('equal')
#         ax.axis('off')
#         plt.show()
    

#     def plot_field(self, u: Callable[[NDArray[Any], NDArray[Any]], NDArray[Any]],
#                    N: int = 100, figsize: Optional[tuple[int, int]] = (16, 2), real_part: bool = False):
#         x = np.linspace(-self.R, self.R, N)
#         y = np.linspace(0., self.H, N)
#         X, Y = np.meshgrid(x, y)
#         Z = u(X, Y)

#         if real_part:
#             Z = np.real(Z)

#         _, ax = plt.subplots(figsize=figsize)

#         ax.pcolorfast((-self.R, self.R), (0., self.H), Z)
#         plt.show()


#     def set_boundary_conditions(self, bc_dict: dict[Region, FluxType]):
#         for region, flux in bc_dict.items():
#             self.domain.edges["flux_type"][self.domain.edges["region"] == region] = flux
#         if np.all(self.domain.edges["flux_type"] >= 0):
#             if self.verbose: 
#                 print('Problem ready for assembly')
#             self.domain.ready_for_assemble = True

class Waveguide(Problem):
    def __init__(self, mesh: TrefftzMesh,
                 boundary_conditions_map: dict[Region, FluxType],
                 k: float, verbose: bool = False,
                 R: float = 5., H: float = 1.):
        self.R = R
        self.H = H
        super().__init__(mesh=mesh, boundary_conditions_map=boundary_conditions_map, verbose=verbose, k=k)

    @classmethod
    def CreateClean(cls, k: float, H: float = 1., R: float = 5., lc: float = 0.3,
                    scatterers: Optional[tuple[Scatterer]] = None, verbose: bool = False):

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

        gmsh.model.addPhysicalGroup(2, [domain], Region.OMEGA, "Omega")
        gmsh.model.addPhysicalGroup(1, [bottom, top], Region.GAMMA, "Gamma")
        gmsh.model.addPhysicalGroup(1, [left, right], Region.SIGMA, "Sigma")
        
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate(2)
        # gmsh.write('CleanWaveguide.msh')
        # self._domain = TrefftzMesh.from_gmsh('CleanWaveguide.msh')
        points, edges, triangles, edges2triangles, locator, cell_sets = GmshArrays(gmsh.model)
        gmsh.finalize()

        mesh = TrefftzMesh(points, edges, triangles, edges2triangles, locator, cell_sets)

        boundary_conditions_map = {
            Region.GAMMA: FluxType.SOUNDHARD,
            Region.SIGMA: FluxType.RADIATING
        }

        return cls(mesh=mesh, boundary_conditions_map=boundary_conditions_map, verbose=verbose, R=R, H=H, k=k)

        # self._regions = Region
    
    # @property
    # def domain(self) -> "TrefftzMesh":
    #     return self._domain
    
    def plot(self, figsize: Optional[tuple[int, int]] = (16, 2), line_width: Optional[int] = 1):
        from matplotlib.collections import LineCollection
        _, ax = plt.subplots(figsize=figsize)
        # ax.triplot(Triangulation(x=M._points[:,0], y=M._points[:,1], triangles=M._triangles),linewidth=lw, color='k')
    
        S = np.where(self.domain.edges["region"] == Region.SIGMA)[0]
        G = np.where(self.domain.edges["region"] == Region.GAMMA)[0]  
        inner = np.where(self.domain.edges["type"] == EdgeType.INNER)[0]
    


        ax.add_collection(LineCollection(np.stack([self.domain.edges[inner]["P"],
                                                   self.domain.edges[inner]["Q"]], axis=1),
                                                   colors='k', linewidths=line_width))

        ax.add_collection(LineCollection(np.stack([self.domain.edges[S]["P"],
                                                   self.domain.edges[S]["Q"]], axis=1),
                                                   colors='r', linewidths=line_width))

        ax.add_collection(LineCollection(np.stack([self.domain.edges[G]["P"],
                                                   self.domain.edges[G]["Q"]], axis=1),
                                                   colors='b', linewidths=line_width))


        ax.axis('equal')
        ax.axis('off')
        plt.show()
    

    def plot_field(self, u: Callable[[NDArray[Any], NDArray[Any]], NDArray[Any]],
                   N: int = 100, figsize: Optional[tuple[int, int]] = (16, 2), real_part: bool = False):
        x = np.linspace(-self.R, self.R, N)
        y = np.linspace(0., self.H, N)
        X, Y = np.meshgrid(x, y)
        Z = u(X, Y)

        if real_part:
            Z = np.real(Z)

        _, ax = plt.subplots(figsize=figsize)

        ax.pcolorfast((-self.R, self.R), (0., self.H), Z)
        ax.axis('equal')
        plt.show()

    def plot_mode(self, n: int):
        self.plot_field(self.mode(n), N=400, real_part=True)

    def mode(self, n: int) -> Callable[[float_array, float_array], complex_array]:
        return Mode(n=n, k=self.k, H=self.H, R=self.R)

    
def CleanWaveguide(k: float, H: float = 1., R: float = 5., lc: float = 0.3,
                scatterers: Optional[tuple[Scatterer]] = None, verbose: bool = False):

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

    gmsh.model.addPhysicalGroup(2, [domain], Region.OMEGA, "Omega")
    gmsh.model.addPhysicalGroup(1, [bottom, top], Region.GAMMA, "Gamma")
    gmsh.model.addPhysicalGroup(1, [left, right], Region.SIGMA, "Sigma")
    
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    # gmsh.write('CleanWaveguide.msh')
    # self._domain = TrefftzMesh.from_gmsh('CleanWaveguide.msh')
    points, edges, triangles, edges2triangles, locator, cell_sets = GmshArrays(gmsh.model)
    gmsh.finalize()

    mesh = TrefftzMesh(points, edges, triangles, edges2triangles, locator, cell_sets)

    boundary_conditions_map = {
        Region.GAMMA: FluxType.SOUNDHARD,
        Region.SIGMA: FluxType.RADIATING
    }

    return Waveguide(mesh=mesh, boundary_conditions_map=boundary_conditions_map, verbose=verbose, R=R, H=H, k=k)



