from trefftz.mesh import TrefftzMesh, BoundaryRegions
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
    def assemble_LHS(self, k: float, edges: NDArray[Any]) -> sparray:
        ...

    def assemble_RHS(self, k: float, edges: NDArray[Any]) -> sparray:
        ...

class SoundHardBC:
    def __init__(self, d_1: float):
        self.d_1 = d_1 # stabilizing parameter
    
    def assemble_LSH(self, k: float, edges: NDArray[Any]) -> sparray:
        pass
        # for edge in edges:
        #     A[] = SoundHard_block(k = k, edge = edge, d = d, d_d = dd, d_1 = self.d_1)

    def assemble_RSH(self, k: float, edges: NDArray[Any]) -> sparray:
        pass
        # for edge in edges:
        #     b[] = SoundHard_block(k = k, edge = edge, d = d, d_d = dd, d_1 = self.d_1)

class RadiatingBC:
    def assemble_LHS(self, k: float, edges: NDArray[Any]) -> sparray:
        ...

    def assemble_RHS(self, k: float, edges: NDArray[Any]) -> sparray:
        ...


class CircularDtN:
    ...

class WaveguideDtN:
    ...
class ExactSolution:
    pass

class BCType(Enum):
    pass
class Problem(Generic[BoundaryRegions]):
    def __init__(self,
                 mesh: TrefftzMesh[BoundaryRegions],
                 wavenumber: float,
                 basis: PlanewaveBasis,
                 boundary_conditions: Mapping[BoundaryRegions, BoundaryCondition],
                 u: ExactSolution | None = None ):
        
        self.mesh = mesh
        self.k = wavenumber
        self.basis = basis
        self.boundary_conditions = boundary_conditions
        self.u = u
        self.A: sparray | None = None
        self.b: sparray | None = None

    @property
    def N_DOF(self):
        return self.basis.N_DOF
    
    @property
    def assembled(self) -> bool:
        return self.A is not None and self.b is not None

    def assemble_LHS(self):
        A = coo_array((self.N_DOF, self.N_DOF), np.complex128)


        # Assembly of the local operators
        data = []
        rows = []
        cols = []
        ## Interior fluxes
        for edge in self.mesh.interior_edges:
            pass 
        
        ## Boundary conditions
        for region in self.mesh.boundary_regions:
            for edge in self.mesh.edges_on_region(region):
                pass


        #     A += self.boundary_conditions[region].assemble_LHS(k = self.k, edges = edges[region], self.N_DOF)
        
        # Assembly of the non-local ones
        A_local = coo_array((data, (rows, cols)), shape=(3, 3))

        
        self.A = A_local + A_nonlocal
    
    def assemble_RHS(self):
        b = coo_array((self.N_DOF), np.complex128)
        for region in self.mesh.boundary_regions:
            print(region)
        #     b += self.boundary_conditions[region].assemble_RHS(k = self.k, edges = edges[region], self.N_DOF)
        self.b = b
    
    def assemble(self):
        self.assemble_RHS()
        self.assemble_LHS()

    def solve(self):
        if self.assembled:
            u_h = spsolve(A=self.A, b=self.b)
        else:
            print("Problem not fully assembled yet.")
        self.u_h = u_h

