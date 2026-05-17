from trefftz.mesh import TrefftzMesh
from trefftz.dg.basis import PlanewaveBasis
from typing import Protocol, TypeVar, Any, Generic
from collections.abc import Mapping
from scipy.sparse.linalg import spsolve
from scipy.sparse import sparray, coo_array
from enum import Enum, EnumMeta
from numpy.typing import NDArray
from trefftz.dg.block_fluxes import SoundHard_block
import numpy as np

class BoundaryCondition(Protocol):
    def assemble_LHS(self, k: float, edges: NDArray[Any], A: sparray) -> None:
        ...

    def assemble_RHS(self, k: float, edges: NDArray[Any], b: sparray) -> None:
        ...

class SoundHardBC:
    def __init__(self, d_1: float):
        self.d_1 = d_1 # stabilizing parameter
    
    def assemble_LSH(self, k: float, edges: NDArray[Any], A: sparray) -> None:
        for edge in edges:
            A[] = SoundHard_block(k = k, edge = edge, d = d, d_d = dd, d_1 = self.d_1)

    def assemble_RSH(self, k: float, edges: NDArray[Any], b: sparray) -> None:
        for edge in edges:
            b[] = SoundHard_block(k = k, edge = edge, d = d, d_d = dd, d_1 = self.d_1)

class RadiatingBC:
    ...


class CircularDtN:
    ...

class WaveguideDtN:
    ...

class BoundaryConditionsConfiguration:
    pass

class ExactSolution:
    pass

BoundaryRegions = TypeVar("BoundaryRegions", bound=Enum)

class Problem(Generic[BoundaryRegions]):
    def __init__(self,
                 mesh: TrefftzMesh[BoundaryRegions],
                 wavenumber: float,
                 basis: PlanewaveBasis,
                 regions: EnumMeta,
                 boundary_conditions: Mapping[BoundaryRegions, BoundaryCondition],
                 u: ExactSolution | None = None ):
        
        self.mesh = mesh
        self.k = wavenumber
        self.basis = basis
        self.regions = regions
        self.boundary_conditions = boundary_conditions
        self.u = u

    @property
    def N_DOF(self):
        return self.basis.N_DOF

    def asseble_LHS(self):
        A = coo_array((self.N_DOF, self.N_DOF), np.complex128)
        for region in self.regions:
            A += self.boundary_conditions[region].assemble_LHS( k = self.k, edges = edges[region], self.N_DOF)
        self.A = A
    
    def asseble_RHS(self):
        b = coo_array((self.N_DOF), np.complex128)
        for region in Regions:
            b += self.boundary_conditions[region].assemble_LHS( k = self.k, edges = edges[region], self.N_DOF)
        self.b = b
    
    def assemble(self):
        self.asseble_RHS()
        self.asseble_LHS()

    def solve(self):
        u_h = spsolve(A=self.A, b=self.b)
        self.u_h = u_h

