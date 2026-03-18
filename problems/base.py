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
    


# # NEW 

# D = TypeVar("D", bound=Domain)
# class Problem(Generic[D]):
#     def __init__(self, domain: D, basis: PlanewaveBasis, k: float, NtD_modes: int, assembler: type(Assemblers),
#                   a: float = 0.5, b: float = 0.5, d_1: float = 0.5, d_2: float = 0.5):
#         self.domain = domain
#         self.basis = basis
#         self.k = k 
#         self.NtD_modes = NtD_modes
#         self.assembler = assembler
#         self.boundary_conditions = None
#         self.stabilizing_parameters = {"a": a, "b": b, "d_1": d_1, "d_2": d_2}

#     def set_boundary_conditions(self, boundary_conditions= Mapping[IntEnum, FluxType]):
#         self.boundary_conditions = boundary_conditions

#     def assembleMatrix(self):
#         if self.boundary_conditions is None:
#             print('no boundary conditions specified, use .set_boundary_conditions')
#             return
#         self.A = SerialAssembleMatrix2(self.domain.mesh.edges, basis=self.basis,
#                                       NtD_modes=self.NtD_modes, boundary_conditions=self.boundary_conditions,
#                                       stabilizing_parameters=self.stabilizing_parameters)

#     def assembleRHS(self, RHS: RSH_type, RHS_params: Mapping[str, int | float]):
#         self.b = SerialAssembleRHS(self.domain.mesh.edges, basis=self.basis,
#                                          NtD_modes=self.NtD_modes, RHS=RHS, RHS_params=RHS_params,
#                                          stabilizing_parameters=self.stabilizing_parameters) ## THE RHS WILL NEED THE BOUNDARY CONDITIONS


#     def assemble(self, RHS: RSH_type, RHS_params: Mapping[str, int | float]):
#         if self.boundary_conditions is None:
#             print('no boundary conditions specified, use .set_boundary_conditions')
#             return
#         self.assembleMatrix()
#         self.assembleRHS(RHS=RHS, RHS_params=RHS_params)
