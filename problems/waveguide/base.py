'''Module for defining a waveguide class, as it will be, I think, the most usefull'''
from trefftz.mesh import TrefftzMesh, EdgeType
# from trefftz.mesh.from_pygmsh import Mesh_from_meshio
from typing import Optional, Callable, Any, Mapping
# from pygmsh.geo import Geometry
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import numpy as np
from numpy.typing import NDArray
# from regions import Region
from problems.base import Domain
import gmsh

from trefftz.mesh.readers import GmshArrays
from problems import Problem
from trefftz.dg.fluxes import FluxType
from problems.waveguide.exact_solutions import WaveguideMode
from trefftz.numpy_types import float_array, complex_array

# from scipy.sparse import csc_matrix
# from trefftz.numpy_types import float_array

from problems.base import AbstractProblem
from enum import IntEnum
from trefftz.dg.basis2 import TrefftzBasis
from trefftz.dg.functions import TrefftzFunction
from scipy.sparse.linalg import spsolve
from .RHS_types import RSH_type

class WaveguideRegions(IntEnum):
    OMEGA = 0
    GAMMA = 1
    SIGMA_L = 3
    SIGMA_R = 4
    SIGMA = 2
    D_OMEGA = 5

# class FluxType(IntEnum):  # should go in flux types or in fluxes
#     TRANSMISSION = 0
#     SOUNDHARD = 1
#     SOUNDSOFT = 2
#     RADIATION = 3

from trefftz.mesh.core import EdgeType
from .assemblers import SerialAssembleMatrix, SerialAssembleRHS, SerialAssembleMatrix2




class WaveguideDomain(Domain[WaveguideRegions]):
    regions = WaveguideRegions

    def __init__(self, mesh: "TrefftzMesh", R: float = 5., H: float = 1.):
        self.R = R
        self.H = H
        self.mesh = mesh

    def plot(self, figsize: tuple[int, int] | None = None, line_width: Optional[int] = 1):

        if figsize is None:
            figsize = 2*int(2*self.R/self.H), 2

        _, ax = plt.subplots(figsize=figsize)
        # ax.triplot(Triangulation(x=M._points[:,0], y=M._points[:,1], triangles=M._triangles),linewidth=lw, color='k')
        regions = self.regions
        edges = self.mesh.edges

        S = np.where(edges["region"] == regions.SIGMA)[0]
        S_L = np.where(edges["region"] == regions.SIGMA_L)[0]
        S_R = np.where(edges["region"] == regions.SIGMA_R)[0]
        G = np.where(edges["region"] == regions.GAMMA)[0]

        inner = np.where(edges["type"] == EdgeType.INNER)[0]  #!!!!!!!!
        ax.add_collection(LineCollection(np.stack([edges[inner]["P"],
                                                   edges[inner]["Q"]], axis=1),
                                                   colors='k', linewidths=line_width))

        colors = {regions.GAMMA: 'b',
                  regions.SIGMA_L: 'r',
                  regions.SIGMA_R: 'g'}

        for reg in colors.keys():
            
            in_reg = np.where(edges["region"] == reg)[0]
            edges_in_reg = edges[in_reg]
            ax.add_collection(LineCollection(np.stack([edges_in_reg["P"], edges_in_reg["Q"]], axis=1),
                                             colors=colors[reg], linewidths=line_width))

        ax.set_title('inner: black, sigma_L: red, sigma_R: green, gamma: blue')
        ax.axis('equal')
        ax.axis('off')
        plt.show()





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

# class Waveguide(Problem):
#     def __init__(self, mesh: TrefftzMesh,
#                  boundary_conditions_map: dict[Region, FluxType],
#                  k: float, verbose: bool = False,
#                  R: float = 5., H: float = 1.):
#         self.R = R
#         self.H = H
#         super().__init__(mesh=mesh, boundary_conditions_map=boundary_conditions_map, verbose=verbose, k=k)

#     @classmethod
#     def CreateClean(cls, k: float, H: float = 1., R: float = 5., lc: float = 0.3,
#                     scatterers: Optional[tuple[Scatterer]] = None, verbose: bool = False):

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
#         gmsh.model.geo.synchronize()

#         gmsh.model.addPhysicalGroup(2, [domain], Region.OMEGA, "Omega")
#         gmsh.model.addPhysicalGroup(1, [bottom, top], Region.GAMMA, "Gamma")
#         gmsh.model.addPhysicalGroup(1, [left, right], Region.SIGMA, "Sigma")
        
#         gmsh.model.geo.synchronize()
#         gmsh.model.mesh.generate(2)
#         # gmsh.write('CleanWaveguide.msh')
#         # self._domain = TrefftzMesh.from_gmsh('CleanWaveguide.msh')
#         points, edges, triangles, edges2triangles, locator, cell_sets = GmshArrays(gmsh.model)
#         gmsh.finalize()

#         mesh = TrefftzMesh(points, edges, triangles, edges2triangles, locator, cell_sets)

#         boundary_conditions_map = {
#             Region.GAMMA: FluxType.SOUNDHARD,
#             Region.SIGMA: FluxType.RADIATING
#         }

#         return cls(mesh=mesh, boundary_conditions_map=boundary_conditions_map, verbose=verbose, R=R, H=H, k=k)

#         # self._regions = Region
    
#     # @property
#     # def domain(self) -> "TrefftzMesh":
#     #     return self._domain
    
#     def plot(self, figsize: Optional[tuple[int, int]] = (16, 2), line_width: Optional[int] = 1):
#         from matplotlib.collections import LineCollection
#         _, ax = plt.subplots(figsize=figsize)
#         # ax.triplot(Triangulation(x=M._points[:,0], y=M._points[:,1], triangles=M._triangles),linewidth=lw, color='k')
    
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
#         ax.axis('equal')
#         plt.show()

#     def plot_mode(self, n: int):
#         self.plot_field(self.mode(n), N=400, real_part=True)

#     def mode(self, n: int) -> Callable[[float_array, float_array], complex_array]:
#         return WaveguideMode(n=n, k=self.k, H=self.H, R=self.R)

    
# def CleanWaveguide(k: float, H: float = 1., R: float = 5., lc: float = 0.3,
#                 scatterers: Optional[tuple[Scatterer]] = None, verbose: bool = False):

#     gmsh.initialize()
#     gmsh.option.setNumber("General.Terminal", int(verbose))
#     gmsh.model.add("Waveguide")
#     p0 = gmsh.model.geo.addPoint(-R, 0., 0., lc)
#     p1 = gmsh.model.geo.addPoint( R, 0., 0., lc)
#     p2 = gmsh.model.geo.addPoint( R,  H, 0., lc)
#     p3 = gmsh.model.geo.addPoint(-R,  H, 0., lc)

#     bottom = gmsh.model.geo.addLine(p0, p1)
#     right  = gmsh.model.geo.addLine(p1, p2)
#     top    = gmsh.model.geo.addLine(p2, p3)
#     left   = gmsh.model.geo.addLine(p3, p0)

#     boundary = gmsh.model.geo.addCurveLoop([bottom, right, top, left])
#     domain = gmsh.model.geo.addPlaneSurface([boundary])
#     gmsh.model.geo.synchronize()

#     gmsh.model.addPhysicalGroup(2, [domain], Region.OMEGA, "Omega")
#     gmsh.model.addPhysicalGroup(1, [bottom, top], Region.GAMMA, "Gamma")
#     gmsh.model.addPhysicalGroup(1, [left, right], Region.SIGMA, "Sigma")
    
#     gmsh.model.geo.synchronize()
#     gmsh.model.mesh.generate(2)
#     # gmsh.write('CleanWaveguide.msh')
#     # self._domain = TrefftzMesh.from_gmsh('CleanWaveguide.msh')
#     points, edges, triangles, edges2triangles, locator, cell_sets = GmshArrays(gmsh.model)
#     gmsh.finalize()

#     mesh = TrefftzMesh(points, edges, triangles, edges2triangles, locator, cell_sets)

#     boundary_conditions_map = {
#         Region.GAMMA: FluxType.SOUNDHARD,
#         Region.SIGMA: FluxType.RADIATING
#     }

#     return Waveguide(mesh=mesh, boundary_conditions_map=boundary_conditions_map, verbose=verbose, R=R, H=H, k=k)


def Clean(H: float = 1., R: float = 5., lc: float = 0.3, verbose: bool = False) -> WaveguideDomain:

    '''Creates a domain corresponging to a waveguide without scatterers.
    It assumes the default tags for the subregions, i.e.:
    - Omega
    - Gamma
    - Sigma
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
    # gmsh.model.addPhysicalGroup(1, [left, right], WaveguideRegions.SIGMA, "Sigma")
    gmsh.model.addPhysicalGroup(1, [left], WaveguideRegions.SIGMA_L, "Sigma_L")
    gmsh.model.addPhysicalGroup(1, [right], WaveguideRegions.SIGMA_R, "Sigma_R")
    
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    # gmsh.write('CleanWaveguide.msh')
    # self._domain = TrefftzMesh.from_gmsh('CleanWaveguide.msh')
    points, edges, triangles, edges2triangles, locator, cell_sets = GmshArrays(gmsh.model)
    gmsh.finalize()

    mesh = TrefftzMesh(points, edges, triangles, edges2triangles, locator, cell_sets)

    return WaveguideDomain(mesh=mesh, R=R, H=H)


# class WaveguideProblem(AbstractProblem[WaveguideDomain, WaveguideRegions]):



#     def plot(self, figsize: Optional[tuple[int, int]] = (16, 2), line_width: Optional[int] = 1):
#         from matplotlib.collections import LineCollection
#         _, ax = plt.subplots(figsize=figsize)
#         # ax.triplot(Triangulation(x=M._points[:,0], y=M._points[:,1], triangles=M._triangles),linewidth=lw, color='k')
    
#         S = np.where(self.domain.mesh.edges["region"] == self.domain.regions.SIGMA)[0]
#         S_L = np.where(self.domain.mesh.edges["region"] == self.domain.regions.SIGMA_L)[0]
#         S_R = np.where(self.domain.mesh.edges["region"] == self.domain.regions.SIGMA_R)[0]
#         G = np.where(self.domain.mesh.edges["region"] == self.domain.regions.GAMMA)[0]
#         # inner = np.where(self.domain.edges["type"] == EdgeType.INNER)[0]  #!!!!!!!!
    


#         # ax.add_collection(LineCollection(np.stack([self.domain.edges[inner]["P"],
#         #                                            self.domain.edges[inner]["Q"]], axis=1),
#         #                                            colors='k', linewidths=line_width))

#         ax.add_collection(LineCollection(np.stack([self.domain.mesh.edges[S]["P"],  
#                                                    self.domain.mesh.edges[S]["Q"]], axis=1),
#                                                    colors='r', linewidths=line_width))

#         ax.add_collection(LineCollection(np.stack([self.domain.mesh.edges[G]["P"],
#                                                    self.domain.mesh.edges[G]["Q"]], axis=1),
#                                                    colors='b', linewidths=line_width))

#         ax.axis('equal')
#         ax.axis('off')
#         plt.show()

#     def plot_field(self, u: Callable[[float_array, float_array], complex_array],
#                    N: int = 100, figsize: Optional[tuple[int, int]] = (16, 2), real_part: bool = False):
#         x = np.linspace(-self.domain.R, self.domain.R, N)
#         y = np.linspace(0., self.domain.H, N)
#         X, Y = np.meshgrid(x, y)
#         Z = u(X, Y)

#         if real_part:
#             Z = np.real(Z)

#         _, ax = plt.subplots(figsize=figsize)

#         ax.pcolorfast((-self.domain.R, self.domain.R), (0., self.domain.H), Z)
#         ax.axis('equal')
#         plt.show()

#     def plot_mode(self, n: int):
#         self.plot_field(self.mode(n), N=400, real_part=True)

#     def mode(self, n: int) -> Callable[[float_array, float_array], complex_array]:
#         return WaveguideMode(n=n, k=self.k, H=self.domain.H, R=self.domain.R)


class Assemblers(IntEnum):
    serial = 1

class WaveguideProblem:
    def __init__(self, domain: WaveguideDomain, N_theta: int, k: float, NtD_modes: int, assembler: type(Assemblers),
                  a: float = 0.5, b: float = 0.5, d_1: float = 0.5, d_2: float = 0.5):
        self.domain = domain
        self.basis = TrefftzBasis(N_elements=domain.mesh.n_triangles, N_theta=N_theta, k=k)
        self.k = k 
        self.NtD_modes = NtD_modes
        self.assembler = assembler
        self.boundary_conditions = None
        self.stabilizing_parameters = {"a": a, "b": b, "d_1": d_1, "d_2": d_2}

    def set_boundary_conditions(self, boundary_conditions= Mapping[WaveguideRegions, FluxType]):
        self.boundary_conditions = boundary_conditions


    def plot_trefftz_function(self, u: TrefftzFunction, figsize: tuple[int, int] | None = None):
        x = np.linspace(-self.domain.R,self.domain.R,200)
        y = np.linspace(0.,self.domain.H, 50)
        X, Y = np.meshgrid(x,y)
        Z = u(X.flatten(), Y.flatten()).reshape(X.shape)
        if figsize is None:
            figsize = 2*int(2*self.domain.R/self.domain.H), 2
        fig, ax = plt.subplots(figsize=figsize)
        ax.pcolormesh(X, Y, np.real(Z), shading='gouraud')
        ax.axis("equal")
        plt.show()

    def plot_field(self, u: Callable[[float_array, float_array], complex_array],
                   N: int = 100, figsize: Optional[tuple[int, int]] | None = None, real_part: bool = False):
        x = np.linspace(-self.domain.R, self.domain.R, N)
        y = np.linspace(0., self.domain.H, N)
        X, Y = np.meshgrid(x, y)
        Z = u(X, Y)

        if figsize is None:
            figsize = 2*int(2*self.domain.R/self.domain.H), 2

        if real_part:
            Z = np.real(Z)

        _, ax = plt.subplots(figsize=figsize)

        ax.pcolorfast((-self.domain.R, self.domain.R), (0., self.domain.H), Z)
        ax.axis('equal')
        plt.show()

    def plot_mode(self, n: int):
        self.plot_field(self.mode(n), N=400, real_part=True)

    def mode(self, n: int) -> Callable[[float_array, float_array], complex_array]:
        return WaveguideMode(n=n, k=self.k, H=self.domain.H, R=self.domain.R)

    def assembleMatrix(self):
        if self.boundary_conditions is None:
            print('no boundary conditions specified, use .set_boundary_conditions')
            return
        self.A = SerialAssembleMatrix2(self.domain.mesh.edges, basis=self.basis,
                                      NtD_modes=self.NtD_modes, boundary_conditions=self.boundary_conditions,
                                      stabilizing_parameters=self.stabilizing_parameters)

    def assembleRHS(self, RHS: RSH_type, RHS_params: Mapping[str, int | float]):
        self.b = SerialAssembleRHS(self.domain.mesh.edges, basis=self.basis,
                                         NtD_modes=self.NtD_modes, RHS=RHS, RHS_params=RHS_params,
                                         stabilizing_parameters=self.stabilizing_parameters) ## THE RHS WILL NEED THE BOUNDARY CONDITIONS


    def assemble(self, RHS: RSH_type, RHS_params: Mapping[str, int | float]):
        if self.boundary_conditions is None:
            print('no boundary conditions specified, use .set_boundary_conditions')
            return
        self.assembleMatrix()
        self.assembleRHS(RHS=RHS, RHS_params=RHS_params)

    def solve(self):
        dofs = spsolve(self.A, self.b)
        u = TrefftzFunction(domain=self.domain, basis=self.basis, dtype=np.complex128)
        u.set(coefficients=dofs)
        return u
         

################# EASY ONE AND THEN BUILD FROM HERE