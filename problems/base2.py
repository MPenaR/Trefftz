from trefftz.mesh import TrefftzMesh, BoundaryRegions
from trefftz.dg.basis import PlanewaveBasis
from typing import Protocol, Generic, Optional, Any
from collections.abc import Mapping
from scipy.sparse.linalg import spsolve
from scipy.sparse import coo_array, csr_array, csc_array
from dataclasses import dataclass
from trefftz.numpy_types import complex_array, float_array, int_array
import numpy as np
from trefftz.dg.serial_kernels import Edge, TT
from trefftz.dg.functions2 import TrefftzFunction

class BoundaryCondition(Protocol):
    data: Optional[Any]


class NeumannBC:
    def __init__(self, data: Any | None = None):
        self.data = data


# class SoundHardBC:
#     def __init__(self, data: object | None = None):
#         self.data = data


class SoundHardBC:
    def __init__(self):
        self.data = None



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
    ...


class Assembler(Protocol):

    def assemble_LHS(self, p: "Problem[BoundaryRegions]") -> coo_array:
        ...

    def assemble_RHS(self, p: "Problem[BoundaryRegions]") -> complex_array:
        ...


class SerialTransmissionKernel(Protocol):
    def LHS(self, edge: Edge, d_phi: float_array, d_psi: float_array, k: float, sign: TT) -> complex:
        ...


class SerialLocalKernel(Protocol):
    def LHS(self, edge: Edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        ...

    def RHS(self, edge: Edge, d_psi: float_array, k: float) -> complex:
        ...


class SerialNonLocalKernel(Protocol):
    def LHS(self, edge_1: Edge, edge_2: Edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        ...


class SerialAssembler(Assembler):

    def assemble_local_bc(self, p: "Problem[BoundaryRegions]", region: BoundaryRegions, rows: list[int], cols: list[int], values: list[complex]):
        bc = p.boundary_conditions[region]
        local_kernel = p.numerics.local_boundary_kernels[type(bc)]
        for edge in p.mesh.edges_on(region):
            T, _ = edge["triangles"]
            for i in p.basis.dofs_on_element(T):
                for j in p.basis.dofs_on_element(T):
                    d_psi = p.basis.global_direction(i)
                    d_phi = p.basis.global_direction(j)
                    value = local_kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, p.k)
                    rows.append(i)
                    cols.append(j)
                    values.append(value)

    def assemble_LHS(self, p: "Problem[BoundaryRegions]") -> coo_array:
        rows: list[int] = []
        cols: list[int] = []
        values: list[complex] = []

        # interior edges (POOR SOLUTION)
        interior_kernel = p.numerics.interior_kernel
        # a = interior_kernel.a 
        # b = interior_kernel.b
        for edge in p.mesh.interior_edges:
            T1, T2 = edge["triangles"]

            # local T1
            # interior_kernel.a = a 
            # interior_kernel.b = b
            sign = TT.PP
            for i in p.basis.dofs_on_element(T1):
                for j in p.basis.dofs_on_element(T1):
                    d_phi = p.basis.global_direction(j)
                    d_psi = p.basis.global_direction(i)
                    val = interior_kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, p.k, sign=sign)
                    rows.append(i)
                    cols.append(j)
                    values.append(val)
            # interior_kernel.a = -a 
            # interior_kernel.b = -b 
            sign = TT.PM
            for i in p.basis.dofs_on_element(T1):
                for j in p.basis.dofs_on_element(T2):
                    d_phi = p.basis.global_direction(j)
                    d_psi = p.basis.global_direction(i)
                    val = interior_kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, p.k, sign=sign)
                    rows.append(i)
                    cols.append(j)
                    values.append(val)

            # interior_kernel.a = a 
            # interior_kernel.b = b
            sign = TT.MP
            for i in p.basis.dofs_on_element(T2):
                for j in p.basis.dofs_on_element(T1):
                    d_phi = p.basis.global_direction(j)
                    d_psi = p.basis.global_direction(i)
                    val = interior_kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, p.k, sign=sign)
                    rows.append(i)
                    cols.append(j)
                    values.append(val)

            # # local T2
            # interior_kernel.a = -a 
            # interior_kernel.b = -b 
            sign = TT.MM
            for i in p.basis.dofs_on_element(T2):
                for j in p.basis.dofs_on_element(T2):
                    d_phi = p.basis.global_direction(j)
                    d_psi = p.basis.global_direction(i)
                    val = interior_kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, p.k, sign=sign)
                    rows.append(i)
                    cols.append(j)
                    values.append(val)

        # boundary conditions implemented as local operators
        for region in p.regions_local_kernel:
            self.assemble_local_bc(p, region, rows, cols, values)


        # boundary conditions implemented as non-local operators
        for region in p.regions_nonlocal_kernel:
            bc = p.boundary_conditions[region]
            non_local_kernel = p.numerics.nonlocal_boundary_kernels[type(bc)]
            for edge_1 in p.mesh.edges_on(region):
                T_1, _ = edge_1["triangles"]
                for edge_2 in p.mesh.edges_on(region):
                    T_2, _ = edge_2["triangles"]
                    for i in p.basis.dofs_on_element(T_1):
                        for j in p.basis.dofs_on_element(T_2):
                            d_phi = p.basis.global_direction(j)
                            d_psi = p.basis.global_direction(i)
                            value = non_local_kernel.LHS(Edge(edge_1["M"], edge_1["l"], edge_1["N"], edge_1["T"]),
                                                         Edge(edge_2["M"], edge_2["l"], edge_2["N"], edge_2["T"]),
                                                         d_phi, d_psi, p.k)
                            rows.append(i)
                            cols.append(j)
                            values.append(value)

        return coo_array((np.asarray(values), (np.asarray(rows), np.asarray(cols))), shape=(p.N_DOF, p.N_DOF))

    def assemble_RHS(self, p: "Problem[BoundaryRegions]") -> complex_array:
        rows: list[int] = []
        values: list[complex] = []
        for region in p.regions_RHS_term:  # I should check redefining this lists as sets or something like that, because of the local AND RHS
            bc = p.boundary_conditions[region]
            local_kernel = p.numerics.local_boundary_kernels[type(bc)]
            # for edge in p.mesh.edges_on_region(region):
            for edge in p.mesh.edges_on(region):
                T, _ = edge["triangles"]
                for i in p.basis.dofs_on_element(T):
                    d_psi = p.basis.global_direction(i)
                    value = local_kernel.RHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_psi, p.k)
                    rows.append(i)
                    values.append(value)

        b = np.zeros((p.N_DOF,), dtype=np.complex128)
        np.add.at(b, rows, values)
        return b


@dataclass
class SerialNumerics:
    interior_kernel: SerialTransmissionKernel
    local_boundary_kernels: Mapping[type[BoundaryCondition], SerialLocalKernel]
    nonlocal_boundary_kernels: Mapping[type[BoundaryCondition], SerialNonLocalKernel]


class Problem(Generic[BoundaryRegions]):
    def __init__(self,
                 mesh: TrefftzMesh[BoundaryRegions],
                 wavenumber: float,
                 basis: PlanewaveBasis,
                 boundary_conditions: Mapping[BoundaryRegions, BoundaryCondition],
                 numerics: SerialNumerics,
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
        self._b: complex_array | None = None

        self._regions_local_kernel = [region for region, bc in self.boundary_conditions.items() if type(bc) in numerics.local_boundary_kernels]
        self._regions_nonlocal_kernel = [region for region, bc in self.boundary_conditions.items() if type(bc) in numerics.nonlocal_boundary_kernels]
        self._regions_RHS_term = [region for region, bc in self.boundary_conditions.items() if bc.data is not None]

    @property
    def regions_local_kernel(self):
        return self._regions_local_kernel

    @property
    def regions_nonlocal_kernel(self):
        return self._regions_nonlocal_kernel
    
    @property
    def regions_RHS_term(self):
        return self._regions_RHS_term

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

    def assemble_RHS(self):
        self._b = self.assembler.assemble_RHS(self)

    @property
    def b(self) -> complex_array | None:
        return self._b


    def assemble(self):
        self.assemble_RHS()
        self.assemble_LHS()

    def solve(self):
        DIRECT = True
        if self.assembled:
            if DIRECT:
                A = csr_array(self.A)
            u_h = TrefftzFunction(self.mesh, self.basis)
            u_h.set(coefficients=spsolve(A=A, b=self.b))
            self.u_h = u_h
        else:
            print("Problem not fully assembled yet.")
