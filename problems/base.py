'''module for defining the problem  base class'''
from typing import TYPE_CHECKING, Optional, Callable, TypeVar, Generic

if TYPE_CHECKING:
    from trefftz.dg.functions import ComplexFunction  # consider changing it to "field"

from trefftz.mesh import TrefftzMesh


from trefftz.numpy_types import float_array, complex_array
from trefftz.dg.fluxes import FluxType
from trefftz.dg.basis import TrefftzBasis, PlanewaveBasis

from enum import IntEnum
from abc import abstractmethod
import numpy as np
from typing import Protocol, Any, Mapping



class Domain:
    mesh: TrefftzMesh
    regions: type[IntEnum]

    def __init__(self, mesh: TrefftzMesh):
        self.mesh = mesh
    
    # @abstractmethod
    # def generate_Trefftzmesh(self):
    #     """Produce TrefftzMesh."""




# class Problem:
#     '''Class which manages defining domains and setting boundary conditions'''

#     def __init__(self, mesh: "TrefftzMesh", boundary_conditions_map: Optional[dict[IntEnum, FluxType]] | None,
#                   k: float, verbose: bool = False):
#         self.domain = mesh
#         self.verbose = verbose
#         if boundary_conditions_map:
#             self.set_boundary_conditions(boundary_conditions_map)
#         self.k = k

#     @abstractmethod
#     def build_geometry(self):
#         """Create gmsh geometry or analytic geometry."""

#     @abstractmethod
#     def generate_Trefftzmesh(self):
#         """Produce TrefftzMesh."""

#     def set_boundary_conditions(self, boundary_conditions_map: dict[IntEnum, FluxType]):
#         self.boundary_conditions_map = boundary_conditions_map
#         """Assign BCs into BoundaryConditionModel."""
#         for region, flux in boundary_conditions_map.items():
#             self.domain.edges["flux_type"][self.domain.edges["region"] == region] = flux
#         if np.all(self.domain.edges["flux_type"] >= 0):
#             if self.verbose: 
#                 print('Problem ready for assembly')
#             self.domain.ready_for_assemble = True

#     @abstractmethod
#     def plot(self, figsize: Optional[tuple[int, int]] = (16, 2), line_width: Optional[int] = 1):
#         """Visualize geometry/mesh/regions."""

#     @abstractmethod
#     def plot_field(self, u: Callable[[float_array, float_array], complex_array],
#                    N: int = 100, figsize: Optional[tuple[int, int]] = (16, 2), real_part: bool = False):
#         """Visualize fields defined on the mesh."""

#     def assemble(self, basis):
#         pass


# class PhysicalModel:
#     def __init__(self, domain: Domain, boundary_conditions, materials):
#         self.domain = domain
#         self.boundary_conditions = boundary_conditions
#         self.materials = materials
#     pass


# class AbstractProblem(Generic[D, R]):
#     '''Class which manages defining domains and setting boundary conditions'''
#     domain: D
#     physics: PhysicalModel
#     verbose: bool
#     regions: R
#     basis: TrefftzBasis
#     k: float


#     def __init__(self, domain: D, model: PhysicalModel, basis: TrefftzBasis, regions: R, k: float, verbose: bool = False):
#         self.domain = domain
#         self.physics = model  # model(domain,k, materials: dict)
#         self.verbose = verbose
#         self.regions = regions
#         self.basis = basis
#         self.k = k

#     @abstractmethod
#     def plot(self, figsize: Optional[tuple[int, int]] = (16, 2), line_width: Optional[int] = 1):
#         """Visualize geometry/mesh/regions."""

#     @abstractmethod
#     def plot_field(self, u: Callable[[float_array, float_array], complex_array],
#                    N: int = 100, figsize: Optional[tuple[int, int]] = (16, 2), real_part: bool = False):
#         """Visualize fields defined on the mesh."""

#     def assemble(self):
#         pass

    # def set_boundary_conditions(self, boundary_conditions_map: dict[IntEnum, FluxType]):
    #     self.boundary_conditions_map = boundary_conditions_map
    #     """Assign BCs into BoundaryConditionModel."""
    #     for region, flux in boundary_conditions_map.items():
    #         self.domain.edges["flux_type"][self.domain.edges["region"] == region] = flux
    #     if np.all(self.domain.edges["flux_type"] >= 0):
    #         if self.verbose:
    #             print('Problem ready for assembly')
    #         self.domain.ready_for_assemble = True



# NEW 

D = TypeVar("D", bound=Domain)
class Problem(Generic[D]):
    def __init__(self, domain: D, basis: PlanewaveBasis, k: float, NtD_modes: int, assembler: type(Assemblers),
                  a: float = 0.5, b: float = 0.5, d_1: float = 0.5, d_2: float = 0.5):
        self.domain = domain
        self.basis = basis
        self.k = k 
        self.NtD_modes = NtD_modes
        self.assembler = assembler
        self.boundary_conditions = None
        self.stabilizing_parameters = {"a": a, "b": b, "d_1": d_1, "d_2": d_2}

    def set_boundary_conditions(self, boundary_conditions= Mapping[IntEnum, FluxType]):
        self.boundary_conditions = boundary_conditions

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
