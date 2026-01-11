'''module for defining the problem  base class'''
from typing import TYPE_CHECKING, Optional, Callable

if TYPE_CHECKING:
    from trefftz.mesh import TrefftzMesh
    from trefftz.dg.functions import ComplexFunction  # consider changing it to "field"

from trefftz.numpy_types import float_array, complex_array
from trefftz.dg.fluxes import FluxType

from enum import IntEnum
from abc import abstractmethod
import numpy as np


class Region(IntEnum):
    pass


class Problem:
    '''Class which manages defining domains and setting boundary conditions'''

    def __init__(self, mesh: "TrefftzMesh", boundary_conditions_map: Optional[dict[Region, FluxType]] | None,
                 verbose: bool = False):
        self.domain = mesh
        self.verbose = verbose
        if boundary_conditions_map:
            self.set_boundary_conditions(boundary_conditions_map)
            



    @abstractmethod
    def build_geometry(self):
        """Create gmsh geometry or analytic geometry."""

    @abstractmethod
    def generate_Trefftzmesh(self):
        """Produce TrefftzMesh."""

    def set_boundary_conditions(self, boundary_conditions_map: dict[Region, FluxType]):
        self.boundary_conditions_map = boundary_conditions_map
        """Assign BCs into BoundaryConditionModel."""
        for region, flux in boundary_conditions_map.items():
            self.domain.edges["flux_type"][self.domain.edges["region"] == region] = flux
        if np.all(self.domain.edges["flux_type"] >= 0):
            if self.verbose: 
                print('Problem ready for assembly')
            self.domain.ready_for_assemble = True

    # def apply_boundary_conditions(self):
    #     """Assign boundary conditions"""
    #     for region, flux in self.boundary_conditions_map.items():
    #         self.domain.edges["flux_type"][self.domain.edges["region"] == region] = flux
    #     if np.all(self.domain.edges["flux_type"] >= 0):
    #         if self.verbose: 
    #             print('Problem ready for assembly')
    #         self.domain.ready_for_assemble = True

    @abstractmethod
    def plot(self, figsize: Optional[tuple[int, int]] = (16, 2), line_width: Optional[int] = 1):
        """Visualize geometry/mesh/regions."""

    @abstractmethod
    def plot_field(self, u: Callable[[float_array, float_array], complex_array],
                   N: int = 100, figsize: Optional[tuple[int, int]] = (16, 2), real_part: bool = False):
        """Visualize fields defined on the mesh."""


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
