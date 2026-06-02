from trefftz.mesh import TrefftzMesh, BoundaryRegions
from trefftz.dg.basis import PlanewaveBasis
from typing import Protocol, TypeVar, Any, Generic
from collections.abc import Mapping
from scipy.sparse.linalg import spsolve
from scipy.sparse import sparray, coo_array, csr_array
from enum import Enum, EnumMeta
from numpy.typing import NDArray
from trefftz.dg.block_fluxes import SoundHard_block
import numpy as np
from trefftz.numpy_types import int_array, complex_array
from numpy.typing import NDArray



class BoundaryCondition(Protocol):
    def assemble_LHS(self, k: float, edge) -> tuple[complex_array, int_array, int_array]:
        ...

    def assemble_RHS(self, k: float, edge) -> tuple[complex_array, int_array]:
        ...

class SoundHardBC:
    def __init__(self, d_1: float):
        self.d_1 = d_1 # stabilizing parameter
    
    def assemble_LSH(self, k: float, edge) -> tuple[complex_array, int_array, int_array]:
        pass
        # for edge in edges:
        #     A[] = SoundHard_block(k = k, edge = edge, d = d, d_d = dd, d_1 = self.d_1)

    def assemble_RSH(self, k: float, edge) -> tuple[complex_array, int_array]:
        pass
        # for edge in edges:
        #     b[] = SoundHard_block(k = k, edge = edge, d = d, d_d = dd, d_1 = self.d_1)

class RadiatingBC:
    def assemble_LHS(self, k: float, edges: NDArray[Any]) -> tuple[np.complex128, int, int]:
        ...

    def assemble_RHS(self, k: float, edges: NDArray[Any]) -> tuple[np.complex128, int]:
        ...


class CircularDtN:
    ...

class WaveguideDtN:
    ...
class ExactSolution:
    pass

class BCType(Enum):
    pass


class Assembler(Protocol):

    def assemble_LSH(self, p: "Problem[BoundaryRegions]") -> sparray:
        ...

    def assemble_RHS(self, p: "Problem[BoundaryRegions]") -> NDArray[Any]
        ...
    
class TransmissionKernel(Protocol):
    def LHS( self, edge, d_phi, d_psi, k_1: float, k_2: float) -> complex:
        ...

class LocalKernel(Protocol):
    def LHS( self, edge, d_phi, d_psi, k: float) -> complex:
        ...
class NonLocalKernel(Protocol):
    def LHS( self, edge_1, edge_2, d_phi, d_psi, k_1: float, k_2: float) -> complex:
        ...


# in reality k_1 and k_2 would belong to phi and psi so it would only require edge1 and edge2

class SerialAssembler(Assembler):

    def assemble_LHS(self, p: "Problem[BoundaryRegions]") -> sparray:
        rows: list[int] = []
        cols: list[int] = []
        data: list[complex] = []
        
        # interior edges
        for edge in p.mesh.interior_edges:
            kernel = p.numerics.interior_kernel
            T1, T2 = edge["triangles"]

            # local T1
            for i in range(p.basis.dofs_on_element(T1)):
                for j in range(p.basis.dofs_on_element(T1)):
                    d_phi = p.basis.D[j]
                    d_psi = p.basis.D[i]
                    val = kernel.LHS(edge, d_phi, d_psi, p.k)
                    rows.append(i)
                    cols.append(j)
                    data.append(val)

            for i in range(p.basis.dofs_on_element(T1)):
                for j in range(p.basis.dofs_on_element(T2)):
                    d_phi = p.basis.D[j]
                    d_psi = p.basis.D[i]
                    val = kernel.LHS(edge, d_phi, d_psi, p.k)
                    rows.append(i)
                    cols.append(j)
                    data.append(val)

            for i in range(p.basis.dofs_on_element(T2)):
                for j in range(p.basis.dofs_on_element(T1)):
                    d_phi = p.basis.D[j]
                    d_psi = p.basis.D[i]
                    val = kernel.LHS(edge, d_phi, d_psi, p.k)
                    rows.append(i)
                    cols.append(j)
                    data.append(val)
            # local T2
            for i in range(p.basis.dofs_on_element(T2)):
                for j in range(p.basis.dofs_on_element(T2)):
                    d_phi = p.basis.D[j]
                    d_psi = p.basis.D[i]
                    val = kernel.LHS(edge, d_phi, d_psi, p.k)
                    rows.append(i)
                    cols.append(j)
                    data.append(val)



        # boundary conditions implemented as local operators
        for region in p.regions_local_kernel:
            bc = p.boundary_conditions[region]
            kernel = p.numerics.local_boundary_kernels[bc]
            for edge in p.mesh.edges_on_region(region):
                T, _ = edge["triangles"]
                for i in p.basis.dofs_on_element(T):
                    for j in p.basis.dofs_on_element(T):
                        d_psi = p.basis.global_direction(i)
                        d_phi = p.basis.global_direction(j)
                        value = kernel.LHS(edge, d_phi, d_psi, p.k)
                        rows.append(i)
                        cols.append(j)
                        data.append(value)

        # boundary conditions implemented as non-local operators
        for region in p.regions_nonlocal_kernel:
            bc = p.boundary_conditions[region]
            kernel = p.numerics.nonlocal_boundary_kernels[bc]
            for edge_1 in p.mesh.edges_on_region(region):
                T_1, _ = edge_1["triangles"]
                for edge_2 in p.mesh.edges_on_region(region):
                    T_2, _ = edge_2["triangles"]
                    for i in range(p.basis.dofs_on_element(T1)):
                        for j in range(p.basis.dofs_on_element(T2)):
                            d_phi = p.basis.D[j]
                            d_psi = p.basis.D[i]
                            row, col, value = kernel.LHS(edge, d_phi, d_psi, p.k)
                            rows.append(row)
                            cols.append(col)
                            data.append(value)
        
        A = coo_array((data, (rows, cols)), shape=(p.N_DOF, p.N_DOF))
        return A.tocsr() #or tocsc depends on iterative or direct
    
    # def assemble_RHS(self, p: "Problem[BoundaryRegions]") -> NDArray[Any]:
    #     rows: list[int] = []
    #     data: list[complex] = []
    #     for edge in p.mesh.boundary_edges:
    #         bc = p.boundary_conditions[edge["region"]]
    #         kernel = SerialKernels[bc]
    #         for i in range(p.basis.N_theta):
    #             d_psi = p.basis.D[i]
    #             row, value = kernel.RHS_data(edge, d_psi, p.k)
    #             rows.append(row)
    #             data.append(value)
    #     b = np.zeros(p.N_DOF, dtype=np.complex128)
    #     np.add.at(b, rows, data)
    #     return b



from dataclasses import dataclass

# @dataclass
# class Numerics:
#     interior_flux: TransmissionKernel
#     local_boundary_fluxes: Mapping[type[BoundaryCondition], LocalKernel]
#     nonlocal_boundary_fluxes: Mapping[type[BoundaryCondition], NonLocalKernel]

@dataclass
class Numerics:
    interior_kernel: TransmissionKernel
    local_boundary_kernels: Mapping[type[BoundaryCondition], LocalKernel]
    nonlocal_boundary_kernels: Mapping[type[BoundaryCondition], NonLocalKernel]

class Problem(Generic[BoundaryRegions]):
    def __init__(self,
                 mesh: TrefftzMesh[BoundaryRegions],
                 wavenumber: float,
                 basis: PlanewaveBasis,
                 boundary_conditions: Mapping[BoundaryRegions, BoundaryCondition],
                 numerics: Numerics,
                 assembler: Assembler = SerialAssembler(),
                 u: ExactSolution | None = None ):
        
        self.mesh = mesh
        self.k = wavenumber
        self.basis = basis
        self.boundary_conditions = boundary_conditions
        self.u = u
        self.numerics = numerics
        self.assembler = assembler
        self.A: sparray | None = None
        self.b: sparray | None = None

        self._regions_local_kernel = [ region for region, bc in self.boundary_conditions.items() if type(bc) in numerics.local_boundary_kernels ]
        self._regions_nonlocal_kernel = [ region for region, bc in self.boundary_conditions.items() if type(bc) in numerics.nonlocal_boundary_kernels ]

    @property
    def regions_local_kernel(self):
        return self._regions_local_kernel

    @property
    def regions_nonlocal_kernel(self):
        return self._regions_nonlocal_kernel


    @property
    def N_DOF(self):
        return self.basis.N_DOF
    
    @property
    def assembled(self) -> bool:
        return self.A is not None and self.b is not None

    # def assemble_LHS(self):

    #     # Assembly of the local operators
    #     data = []
    #     rows = []
    #     cols = []
    #     ## Interior fluxes
    #     for edge in self.mesh.interior_edges:
    #         data_edge, rows_edge, cols_edge = TransmissionFlux.assemble(edge)
    #         data.extend(data_edge)
    #         rows.extend(rows_edge)
    #         cols.extend(cols_edge)
        
    #     ## Boundary conditions
    #     for region in self.mesh.boundary_regions:
    #         for edge in self.mesh.edges_on_region(region):
    #                         data_edge, rows_edge, cols_edge = self.boundary_conditions[region].assemble_LHS(self.k, edge)
    #         data.extend(data_edge)
    #         rows.extend(rows_edge)
    #         cols.extend(cols_edge)

    #     # Assembly of the local ones
    #     A_local = coo_array((data, (rows, cols)), shape=(self.N_DOF, self.N_DOF))

    #     #     A += self.boundary_conditions[region].assemble_LHS(k = self.k, edges = edges[region], self.N_DOF)
        
    #     A_non_local = coo_array((self.N_DOF, self.N_DOF), dtype=np.complex128)
    #     self.A: csr_array = (A_local + A_non_local).tocsr()
    
    # def assemble_RHS(self):
    #     b = coo_array((self.N_DOF), np.complex128)
    #     for region in self.mesh.boundary_regions:
    #         print(region)
    #     #     b += self.boundary_conditions[region].assemble_RHS(k = self.k, edges = edges[region], self.N_DOF)
    #     self.b = b

    def assemble_RHS(self):
        self.A = self.assembler.assemble_RHS(self)

    def assemble_LHS(self):
        self.b = self.assembler.assemble_LHS(self)

    def assemble(self):
        self.assemble_RHS()
        self.assemble_LHS()

    def solve(self):
        if self.assembled:
            u_h = spsolve(A=self.A, b=self.b)
        else:
            print("Problem not fully assembled yet.")
        self.u_h = u_h


class SoundHardKernel:
    '''Serial SoundHard kernel'''
    def __init__(self, d_1: float):
        self.d_1 = d_1
    
    def LSH(self, edge, d_phi, d_psi, k: float) -> complex:
        d_1 = self.d_1
        d_m = d_psi
        d_n = d_phi

        M = edge.M
        l = edge.l
        N = edge.N
        T = edge.T

        return -1j*k*l*(1 + d_1 * dot(d_n, N))*dot(d_m, N)*exp(1j*k*dot(d_n - d_m, M)) * sinc(k*l/(2*pi)*dot(d_n-d_m, T))
    
class UltraWeakKernel:
    '''Transmission kernel fro the UWVF'''
    def __init__(self, a: float, b: float):
        self.a = a 
        self.b = b
    
    def LSH(self, edge, d_phi, d_psi, k_1: float, k_2: float) -> complex:
        a = self.a 
        b = self.b
        d_m = d_psi
        d_n = d_phi

        k_n = k_1 
        k_m = k_2


        M = edge.M
        N = edge.N
        T = edge.T
        l = edge.l

        # I = -1j*l/2*(2*a*k + k_n*dot(d_n, N) + k_m*dot(d_m, N) + 2*b/k*k_n*dot(d_n, N)*k_m*dot(d_m, N))*exp(1j*dot(k_n*d_n - k_m*d_m, M))*sinc(l/(2*pi)*dot(k_n*d_n - k_m*d_m,T))
        return -1j*k*l*( 1/2*( k_n/k*dot(d_n, N) + k_m/k*dot(d_m, N)) + a + b*k_n/k*k_m/k*dot(d_n, N)*dot(d_m, N))*exp(1j*dot(k_n*d_n - k_m*d_m, M))*sinc(l/(2*pi)*dot(k_n*d_n - k_m*d_m,T))
    




