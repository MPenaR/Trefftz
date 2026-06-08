from trefftz.mesh import TrefftzMesh, BoundaryRegions
from trefftz.dg.basis import PlanewaveBasis
from typing import Protocol, Generic, Optional
from collections.abc import Mapping
from scipy.sparse.linalg import spsolve
from scipy.sparse import coo_array# , csr_array, csc_array
from dataclasses import dataclass
from trefftz.numpy_types import complex_array
import numpy as np


class BoundaryCondition(Protocol):
    data: Optional[object]


class SoundHardBC:
    def __init__(self, data: object | None = None):
        self.data = data

class RadiatingBC:
    def __init__(self, data: object | None = None):
        self.data = data

class NtDBC:
    def __init__(self, truncating_radius: float, data: object | None = None):
        self.truncating_radius = truncating_radius
        self.data = data
    


class CircularNtD:
    ...

class WaveguideNtD:
    ...
class ExactSolution:
    pass


class Assembler(Protocol):

    def assemble_LHS(self, p: "Problem[BoundaryRegions]") -> coo_array:
        ...

    # def assemble_RHS(self, p: "Problem[BoundaryRegions]") -> NDArray[Any]
    #     ...
    

# not a good protocol right now, it should expect k, eps_1 and eps_2, or neither of them as they belong to the basis    
class TransmissionKernel(Protocol):
    def LHS( self, edge, d_phi, d_psi, k: float) -> complex:
        ...

class LocalKernel(Protocol):
    def LHS( self, edge, d_phi, d_psi, k: float) -> complex:
        ...
class NonLocalKernel(Protocol):
    def LHS( self, edge_1, edge_2, d_phi, d_psi, k: float) -> complex:
        ...

class SerialAssembler(Assembler):

    def assemble_LHS(self, p: "Problem[BoundaryRegions]") -> coo_array:
        rows: list[int] = []
        cols: list[int] = []
        values: list[complex] = []
        
        # interior edges
        for edge in p.mesh.interior_edges:
            interior_kernel = p.numerics.interior_kernel
            T1, T2 = edge["triangles"]

            # local T1
            for i in p.basis.dofs_on_element(T1):
                for j in p.basis.dofs_on_element(T1):
                    d_phi = p.basis.global_direction(j)
                    d_psi = p.basis.global_direction(i)
                    val = interior_kernel.LHS(edge, d_phi, d_psi, p.k)
                    rows.append(i)
                    cols.append(j)
                    values.append(val)

            for i in p.basis.dofs_on_element(T1):
                for j in p.basis.dofs_on_element(T2):
                    d_phi = p.basis.global_direction(j)
                    d_psi = p.basis.global_direction(i)
                    val = interior_kernel.LHS(edge, d_phi, d_psi, p.k)
                    rows.append(i)
                    cols.append(j)
                    values.append(val)

            for i in p.basis.dofs_on_element(T2):
                for j in p.basis.dofs_on_element(T1):
                    d_phi = p.basis.global_direction(j)
                    d_psi = p.basis.global_direction(i)
                    val = interior_kernel.LHS(edge, d_phi, d_psi, p.k)
                    rows.append(i)
                    cols.append(j)
                    values.append(val)
            # local T2
            for i in p.basis.dofs_on_element(T2):
                for j in p.basis.dofs_on_element(T2):
                    d_phi = p.basis.global_direction(j)
                    d_psi = p.basis.global_direction(i)
                    val = interior_kernel.LHS(edge, d_phi, d_psi, p.k)
                    rows.append(i)
                    cols.append(j)
                    values.append(val)



        # boundary conditions implemented as local operators
        for region in p.regions_local_kernel:
            bc = p.boundary_conditions[region]
            local_kernel = p.numerics.local_boundary_kernels[type(bc)]
            for edge in p.mesh.edges_on_region(region):
                T, _ = edge["triangles"]
                for i in p.basis.dofs_on_element(T):
                    for j in p.basis.dofs_on_element(T):
                        d_psi = p.basis.global_direction(i)
                        d_phi = p.basis.global_direction(j)
                        value = local_kernel.LHS(edge, d_phi, d_psi, p.k)
                        rows.append(i)
                        cols.append(j)
                        values.append(value)

        # boundary conditions implemented as non-local operators
        for region in p.regions_nonlocal_kernel:
            bc = p.boundary_conditions[region]
            non_local_kernel = p.numerics.nonlocal_boundary_kernels[type(bc)]
            for edge_1 in p.mesh.edges_on_region(region):
                T_1, _ = edge_1["triangles"]
                for edge_2 in p.mesh.edges_on_region(region):
                    T_2, _ = edge_2["triangles"]
                    for i in p.basis.dofs_on_element(T_1):
                        for j in p.basis.dofs_on_element(T_2):
                            d_phi = p.basis.global_direction(j)
                            d_psi = p.basis.global_direction(i)
                            value = non_local_kernel.LHS(edge_1, edge_2, d_phi, d_psi, p.k)
                            rows.append(i)
                            cols.append(j)
                            values.append(value)

        return coo_array((values, (rows, cols)), shape=(p.N_DOF, p.N_DOF))
    
    def assemble_RHS(self, p: "Problem[BoundaryRegions]") -> complex_array:
        rows: list[int] = []
        values: list[complex] = []
        
        # boundary conditions implemented as local operators
        for region in p.regions_local_kernel:
            bc = p.boundary_conditions[region]
            if bc.data is None:
                continue
            local_kernel = p.numerics.local_boundary_kernels[type(bc)]
            for edge in p.mesh.edges_on_region(region):
                T, _ = edge["triangles"]
                for i in p.basis.dofs_on_element(T):
                    d_psi = p.basis.global_direction(i)
                    value = local_kernel.RHS(edge, d_psi, p.k)
                    rows.append(i)
                    values.append(value)

        # boundary conditions implemented as non-local operators
        for region in p.regions_nonlocal_kernel:
            bc = p.boundary_conditions[region]
            if bc.data is None:
                continue
            non_local_kernel = p.numerics.nonlocal_boundary_kernels[type(bc)]
            for edge_1 in p.mesh.edges_on_region(region):
                T_1, _ = edge_1["triangles"]
                for i in p.basis.dofs_on_element(T_1):
                    d_psi = p.basis.global_direction(i)
                    value = non_local_kernel.RHS(edge_1, d_psi, p.k)
                    rows.append(i)
                    values.append(value)
        b = np.zeros((p.N_DOF,), dtype=np.complex128)
        np.add.at(b, rows, values)
        return b
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
        self._A: coo_array | None = None
        self._b: complex_array| None = None

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

    def assemble_LHS(self):
        self._A = self.assembler.assemble_LHS(self)
        return self._A

    @property
    def A(self) -> coo_array | None:
        return self._A

    @property
    def b(self) -> complex_array | None:
        return self._b


    def assemble_RHS(self):
        pass
        # self.A = self.assembler.assemble_RHS(self)

    def assemble(self):
        self.assemble_RHS()
        self.assemble_LHS()

    def solve(self):
        if self.assembled:
            u_h = spsolve(A=self.A, b=self.b)
            self.u_h = u_h
        else:
            print("Problem not fully assembled yet.")

    




