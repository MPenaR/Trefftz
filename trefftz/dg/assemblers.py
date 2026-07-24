from trefftz.mesh import TrefftzMesh
from scipy.sparse import coo_array, bsr_array, csr_array, csc_array
from trefftz.dg.serial_kernels import Edge, SIGN
import numpy as np
from trefftz.numpy_types import complex_array
from trefftz.dg.serial_kernels import SerialLocalKernel, SerialNonLocalKernel, SerialTransmissionKernel
from trefftz.dg.block_kernels import BlockLocalKernel, BlockNonLocalKernel, BlockTransmissionKernel
from trefftz.dg.basis import PlanewaveBasis
from typing import Protocol, Mapping, Any
from dataclasses import dataclass
from trefftz.dg.boundary_conditions import BoundaryCondition
from abc import ABC, abstractmethod
from enum import IntEnum, StrEnum

## THIS IS WRONG, CHECK IT LATER FOR THE BLOCK KERNELS 
@dataclass
class Numerics(Protocol):
    interior_kernel: SerialTransmissionKernel
    local_boundary_kernels: Mapping[type[BoundaryCondition], SerialLocalKernel]
    nonlocal_boundary_kernels: Mapping[type[BoundaryCondition], SerialNonLocalKernel]


@dataclass
class SerialNumerics:
    interior_kernel: SerialTransmissionKernel
    local_boundary_kernels: Mapping[type[BoundaryCondition], SerialLocalKernel]
    nonlocal_boundary_kernels: Mapping[type[BoundaryCondition], SerialNonLocalKernel]


@dataclass
class BlockNumerics:
    interior_kernel: BlockTransmissionKernel
    local_boundary_kernels: Mapping[type[BoundaryCondition], BlockLocalKernel]
    nonlocal_boundary_kernels: Mapping[type[BoundaryCondition], BlockNonLocalKernel]


class Assembler[BR: (StrEnum, IntEnum), num: Numerics](ABC):
    def __init__(self,
                 mesh: TrefftzMesh[BR],
                 boundary_conditions: Mapping[BR, BoundaryCondition],
                 numerics: num,
                 basis: PlanewaveBasis):
        
        self._mesh = mesh
        self._boundary_conditions = boundary_conditions
        self._numerics = numerics
        self._basis = basis
        self._regions_local_kernel = [region for region, bc in boundary_conditions.items() if type(bc) in numerics.local_boundary_kernels]
        self._regions_nonlocal_kernel = [region for region, bc in boundary_conditions.items() if type(bc) in numerics.nonlocal_boundary_kernels]
        self._regions_RHS_term = [region for region, bc in boundary_conditions.items() if bc.data is not None]


    @abstractmethod
    def assemble_LHS(self) -> coo_array:
        ...

    @abstractmethod
    def assemble_RHS(self) -> complex_array:
        ...


class SerialAssembler(Assembler[Any, SerialNumerics]):

    def assemble_local_bc(self,
                          edges_on_region,
                          kernel: SerialLocalKernel,
                          basis: PlanewaveBasis,
                          rows: list[int],
                          cols: list[int],
                          values: list[complex]):
   
        for edge in edges_on_region:
            T, _ = edge["triangles"]
            for i in basis.dofs_on_element(T):
                for j in basis.dofs_on_element(T):
                    d_psi = basis.global_direction(i)
                    d_phi = basis.global_direction(j)
                    value = kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, basis.k)
                    rows.append(i)
                    cols.append(j)
                    values.append(value)

    def assemble_LHS(self) -> coo_array:

        regions_local_kernel = self._regions_local_kernel
        regions_nonlocal_kernel = self._regions_nonlocal_kernel
        mesh = self._mesh
        boundary_conditions = self._boundary_conditions
        numerics = self._numerics
        basis = self._basis


        rows: list[int] = []
        cols: list[int] = []
        values: list[complex] = []

        # interior edges
        interior_kernel = numerics.interior_kernel
        for edge in mesh.interior_edges:
            for (i_v, T_v) in enumerate(edge["triangles"]):
                for (i_u, T_u) in enumerate(edge["triangles"]):
                    sign = SIGN((i_u, i_v))
                    for i in basis.dofs_on_element(T_v):
                        for j in basis.dofs_on_element(T_u):
                            d_phi = basis.global_direction(j)
                            d_psi = basis.global_direction(i)
                            val = interior_kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, basis.k, sign=sign)
                            rows.append(i)
                            cols.append(j)
                            values.append(val)

        # boundary conditions implemented as local operators
        for region in regions_local_kernel:
            edges_on_region = mesh.edges_on(region)
            bc = boundary_conditions[region]
            kernel = numerics.local_boundary_kernels[type(bc)]
            self.assemble_local_bc(edges_on_region, kernel, basis, rows, cols, values)


        # boundary conditions implemented as non-local operators
        for region in regions_nonlocal_kernel:
            bc = boundary_conditions[region]
            non_local_kernel = numerics.nonlocal_boundary_kernels[type(bc)]
            for edge_1 in mesh.edges_on(region):
                T_1, _ = edge_1["triangles"]
                for edge_2 in mesh.edges_on(region):
                    T_2, _ = edge_2["triangles"]
                    for i in basis.dofs_on_element(T_1):
                        for j in basis.dofs_on_element(T_2):
                            d_phi = basis.global_direction(j)
                            d_psi = basis.global_direction(i)
                            value = non_local_kernel.LHS(Edge(edge_1["M"], edge_1["l"], edge_1["N"], edge_1["T"]),
                                                         Edge(edge_2["M"], edge_2["l"], edge_2["N"], edge_2["T"]),
                                                         d_phi, d_psi, basis.k)
                            rows.append(i)
                            cols.append(j)
                            values.append(value)

        return coo_array((np.asarray(values), (np.asarray(rows), np.asarray(cols))), shape=(basis.N_DOF, basis.N_DOF))

    def assemble_RHS(self) -> complex_array:

        rows: list[int] = []
        values: list[complex] = []

        regions_RHS_term = self._regions_RHS_term
        mesh = self._mesh
        boundary_conditions = self._boundary_conditions
        numerics = self._numerics
        basis = self._basis

        for region in regions_RHS_term:  # I should check redefining this lists as sets or something like that, because of the local AND RHS
            bc = boundary_conditions[region]
            local_kernel = numerics.local_boundary_kernels[type(bc)]
            # for edge in p.mesh.edges_on_region(region):
            for edge in mesh.edges_on(region):
                T, _ = edge["triangles"]
                for i in basis.dofs_on_element(T):
                    d_psi = basis.global_direction(i)
                    value = local_kernel.RHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_psi, basis.k)
                    rows.append(i)
                    values.append(value)

        b = np.zeros((basis.N_DOF,), dtype=np.complex128)
        np.add.at(b, rows, values)
        return b


class BlockAssembler(Assembler[Any, BlockNumerics]):

    # def assemble_local_bc(self,
    #                       edges_on_region,
    #                       kernel: SerialLocalKernel,
    #                       basis: TrefftzBasis,
    #                       rows: list[int],
    #                       cols: list[int],
    #                       values: list[complex]):
   
    #     for edge in edges_on_region:
    #         T, _ = edge["triangles"]
    #         for i in basis.dofs_on_element(T):
    #             for j in basis.dofs_on_element(T):
    #                 d_psi = basis.global_direction(i)
    #                 d_phi = basis.global_direction(j)
    #                 value = kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, basis.k)
    #                 rows.append(i)
    #                 cols.append(j)
    #                 values.append(value)

    def assemble_LHS(self) -> coo_array:
        raise NotImplementedError("Yet")

        regions_local_kernel = self._regions_local_kernel
        regions_nonlocal_kernel = self._regions_nonlocal_kernel
        mesh = self._mesh
        boundary_conditions = self._boundary_conditions
        numerics = self._numerics
        basis = self._basis


        # interior edges
        interior_kernel = numerics.interior_kernel
        for edge in mesh.interior_edges:
            for (i_v, T_v) in enumerate(edge["triangles"]):
                for (i_u, T_u) in enumerate(edge["triangles"]):
                    sign = SIGN((i_u, i_v))
                    for i in basis.dofs_on_element(T_v):
                        for j in basis.dofs_on_element(T_u):
                            d_phi = basis.global_direction(j)
                            d_psi = basis.global_direction(i)
                            val = interior_kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, basis.k, sign=sign)
                            rows.append(i)
                            cols.append(j)
                            values.append(val)


        # boundary conditions implemented as local operators  #this should be generalized for a plane wave basis with non homogeneous number of waves
        # that imples that basis.D is a function on the element ID as well as D_D
        # only for homogeneous P makes sense to precompute D_D as it is reused everywhere. For non homogeneous D_D there is a D_D per pair of T's and
        # they are not reused in general.
        for region in regions_local_kernel:
            bc = boundary_conditions[region]
            kernel = numerics.local_boundary_kernels[type(bc)]
            for edge in mesh.edges_on(region):
                T, _ = edge["triangles"]
                block = kernel.LHS(edge, basis.D, basis.D_D, basis.k)
                


        # boundary conditions implemented as non-local operators
        for region in regions_nonlocal_kernel:
            bc = boundary_conditions[region]
            non_local_kernel = numerics.nonlocal_boundary_kernels[type(bc)]
            for edge_1 in mesh.edges_on(region):
                T_1, _ = edge_1["triangles"]
                for edge_2 in mesh.edges_on(region):
                    T_2, _ = edge_2["triangles"]
                    for i in basis.dofs_on_element(T_1):
                        for j in basis.dofs_on_element(T_2):
                            d_phi = basis.global_direction(j)
                            d_psi = basis.global_direction(i)
                            value = non_local_kernel.LHS(Edge(edge_1["M"], edge_1["l"], edge_1["N"], edge_1["T"]),
                                                         Edge(edge_2["M"], edge_2["l"], edge_2["N"], edge_2["T"]),
                                                         d_phi, d_psi, basis.k)
                            rows.append(i)
                            cols.append(j)
                            values.append(value)

        return coo_array((np.asarray(values), (np.asarray(rows), np.asarray(cols))), shape=(basis.N_DOF, basis.N_DOF))

    def assemble_RHS(self) -> complex_array:

        rows: list[int] = []
        values: list[complex] = []

        regions_RHS_term = self._regions_RHS_term
        mesh = self._mesh
        boundary_conditions = self._boundary_conditions
        numerics = self._numerics
        basis = self._basis

        for region in regions_RHS_term:  # I should check redefining this lists as sets or something like that, because of the local AND RHS
            bc = boundary_conditions[region]
            local_kernel = numerics.local_boundary_kernels[type(bc)]
            # for edge in p.mesh.edges_on_region(region):
            for edge in mesh.edges_on(region):
                T, _ = edge["triangles"]
                for i in basis.dofs_on_element(T):
                    d_psi = basis.global_direction(i)
                    value = local_kernel.RHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_psi, basis.k)
                    rows.append(i)
                    values.append(value)

        b = np.zeros((basis.N_DOF,), dtype=np.complex128)
        np.add.at(b, rows, values)
        return b
