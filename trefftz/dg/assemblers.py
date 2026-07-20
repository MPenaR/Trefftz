from trefftz.mesh import TrefftzMesh, BoundaryRegions
from scipy.sparse import coo_array, csr_array, csc_array
from trefftz.dg.serial_kernels import Edge, TT
import numpy as np
from trefftz.numpy_types import complex_array
from trefftz.dg.serial_kernels import SerialLocalKernel, SerialNonLocalKernel, SerialTransmissionKernel
from trefftz.dg.block_kernels import BlockLocalKernel
from trefftz.dg.basis import TrefftzBasis
from typing import Protocol


class Assembler(Protocol):

    def assemble_LHS(self, p: "Problem[BoundaryRegions]", k: float) -> coo_array:
        ...

    def assemble_RHS(self, p: "Problem[BoundaryRegions]", k: float) -> complex_array:
        ...


class SerialAssembler(Assembler):

    def assemble_local_bc(self,
                          edges_on_region,
                          kernel: SerialLocalKernel,
                          basis: TrefftzBasis,
                          rows: list[int],
                          cols: list[int],
                          values: list[complex]):
   
        for edge in edges_on_region:
            T, _ = edge["triangles"]
            for i in basis.dofs_on_element(T):
                for j in basis.dofs_on_element(T):
                    d_psi = basis.global_direction(i)
                    d_phi = basis.global_direction(j)
                    value = kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, p.k)
                    rows.append(i)
                    cols.append(j)
                    values.append(value)

    # def assemble_LHS(self, p: "Problem[BoundaryRegions]") -> coo_array:
    def assemble_LHS(self,
                     numerics,
                     mesh: TrefftzMesh[BoundaryRegions],
                     basis: TrefftzBasis,
                     ) -> coo_array:

        rows: list[int] = []
        cols: list[int] = []
        values: list[complex] = []

        # interior edges (POOR SOLUTION)
        interior_kernel = numerics.interior_kernel
        # a = interior_kernel.a 
        # b = interior_kernel.b
        for edge in mesh.interior_edges:
            T1, T2 = edge["triangles"]

            # local T1
            # interior_kernel.a = a 
            # interior_kernel.b = b
            sign = TT.PP
            for i in basis.dofs_on_element(T1):
                for j in basis.dofs_on_element(T1):
                    d_phi = basis.global_direction(j)
                    d_psi = basis.global_direction(i)
                    val = interior_kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, p.k, sign=sign)
                    rows.append(i)
                    cols.append(j)
                    values.append(val)
            # interior_kernel.a = -a 
            # interior_kernel.b = -b 
            sign = TT.PM
            for i in basis.dofs_on_element(T1):
                for j in basis.dofs_on_element(T2):
                    d_phi = basis.global_direction(j)
                    d_psi = basis.global_direction(i)
                    val = interior_kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, p.k, sign=sign)
                    rows.append(i)
                    cols.append(j)
                    values.append(val)

            # interior_kernel.a = a 
            # interior_kernel.b = b
            sign = TT.MP
            for i in basis.dofs_on_element(T2):
                for j in basis.dofs_on_element(T1):
                    d_phi = basis.global_direction(j)
                    d_psi = basis.global_direction(i)
                    val = interior_kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, p.k, sign=sign)
                    rows.append(i)
                    cols.append(j)
                    values.append(val)

            # # local T2
            # interior_kernel.a = -a 
            # interior_kernel.b = -b 
            sign = TT.MM
            for i in basis.dofs_on_element(T2):
                for j in basis.dofs_on_element(T2):
                    d_phi = basis.global_direction(j)
                    d_psi = basis.global_direction(i)
                    val = interior_kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, p.k, sign=sign)
                    rows.append(i)
                    cols.append(j)
                    values.append(val)

        # boundary conditions implemented as local operators
        for region in p.regions_local_kernel:
            edges_on_region = mesh.edges_on(region)
            bc = p.boundary_conditions[region]
            kernel = numerics.local_boundary_kernels[type(bc)]
            self.assemble_local_bc(edges_on_region, kernel, basis, rows, cols, values)


        # boundary conditions implemented as non-local operators
        for region in p.regions_nonlocal_kernel:
            bc = p.boundary_conditions[region]
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
                                                         d_phi, d_psi, p.k)
                            rows.append(i)
                            cols.append(j)
                            values.append(value)

        return coo_array((np.asarray(values), (np.asarray(rows), np.asarray(cols))), shape=(p.N_DOF, p.N_DOF))

    # def assemble_RHS(self, p: "Problem[BoundaryRegions]") -> complex_array:
    def assemble_RHS(self,
                     mesh: TrefftzMesh[BoundaryRegions],
                     boundary_conditions,
                     numerics,
                     basis: TrefftzBasis) -> complex_array:

        rows: list[int] = []
        values: list[complex] = []
        for region in p.regions_RHS_term:  # I should check redefining this lists as sets or something like that, because of the local AND RHS
            bc = boundary_conditions[region]
            local_kernel = numerics.local_boundary_kernels[type(bc)]
            # for edge in p.mesh.edges_on_region(region):
            for edge in mesh.edges_on(region):
                T, _ = edge["triangles"]
                for i in basis.dofs_on_element(T):
                    d_psi = basis.global_direction(i)
                    value = local_kernel.RHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_psi, p.k)
                    rows.append(i)
                    values.append(value)

        b = np.zeros((basis.N_DOF,), dtype=np.complex128)
        np.add.at(b, rows, values)
        return b



class BlockAssembler(Assembler):

    def assemble_local_bc(self,
                          edges_on_region,
                          kernel: BlockLocalKernel,
                          basis: TrefftzBasis
                          k: float):
   
        for edge in edges_on_region:
            T, _ = edge["triangles"]
            block = kernel.LHS(edge, basis.d, basis.D, k)
            for i in basis.dofs_on_element(T):
                for j in basis.dofs_on_element(T):
                    d_psi = basis.global_direction(i)
                    d_phi = basis.global_direction(j)
                    value = kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, p.k)
                    rows.append(i)
                    cols.append(j)
                    values.append(value)

    # def assemble_LHS(self, p: "Problem[BoundaryRegions]") -> coo_array:
    def assemble_LHS(self,
                     numerics,
                     mesh: TrefftzMesh[BoundaryRegions],
                     basis: TrefftzBasis,
                     ) -> coo_array:

        rows: list[int] = []
        cols: list[int] = []
        values: list[complex] = []

        # interior edges (POOR SOLUTION)
        interior_kernel = numerics.interior_kernel
        # a = interior_kernel.a 
        # b = interior_kernel.b
        for edge in mesh.interior_edges:
            T1, T2 = edge["triangles"]

            # local T1
            # interior_kernel.a = a 
            # interior_kernel.b = b
            sign = TT.PP
            for i in basis.dofs_on_element(T1):
                for j in basis.dofs_on_element(T1):
                    d_phi = basis.global_direction(j)
                    d_psi = basis.global_direction(i)
                    val = interior_kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, p.k, sign=sign)
                    rows.append(i)
                    cols.append(j)
                    values.append(val)
            # interior_kernel.a = -a 
            # interior_kernel.b = -b 
            sign = TT.PM
            for i in basis.dofs_on_element(T1):
                for j in basis.dofs_on_element(T2):
                    d_phi = basis.global_direction(j)
                    d_psi = basis.global_direction(i)
                    val = interior_kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, p.k, sign=sign)
                    rows.append(i)
                    cols.append(j)
                    values.append(val)

            # interior_kernel.a = a 
            # interior_kernel.b = b
            sign = TT.MP
            for i in basis.dofs_on_element(T2):
                for j in basis.dofs_on_element(T1):
                    d_phi = basis.global_direction(j)
                    d_psi = basis.global_direction(i)
                    val = interior_kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, p.k, sign=sign)
                    rows.append(i)
                    cols.append(j)
                    values.append(val)

            # # local T2
            # interior_kernel.a = -a 
            # interior_kernel.b = -b 
            sign = TT.MM
            for i in basis.dofs_on_element(T2):
                for j in basis.dofs_on_element(T2):
                    d_phi = basis.global_direction(j)
                    d_psi = basis.global_direction(i)
                    val = interior_kernel.LHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_phi, d_psi, p.k, sign=sign)
                    rows.append(i)
                    cols.append(j)
                    values.append(val)

        # boundary conditions implemented as local operators
        for region in p.regions_local_kernel:
            edges_on_region = mesh.edges_on(region)
            bc = p.boundary_conditions[region]
            kernel = numerics.local_boundary_kernels[type(bc)]
            self.assemble_local_bc(edges_on_region, kernel, basis, rows, cols, values)


        # boundary conditions implemented as non-local operators
        for region in p.regions_nonlocal_kernel:
            bc = p.boundary_conditions[region]
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
                                                         d_phi, d_psi, p.k)
                            rows.append(i)
                            cols.append(j)
                            values.append(value)

        return coo_array((np.asarray(values), (np.asarray(rows), np.asarray(cols))), shape=(p.N_DOF, p.N_DOF))

    # def assemble_RHS(self, p: "Problem[BoundaryRegions]") -> complex_array:
    def assemble_RHS(self,
                     mesh: TrefftzMesh[BoundaryRegions],
                     boundary_conditions,
                     numerics,
                     basis: TrefftzBasis) -> complex_array:

        rows: list[int] = []
        values: list[complex] = []
        for region in p.regions_RHS_term:  # I should check redefining this lists as sets or something like that, because of the local AND RHS
            bc = boundary_conditions[region]
            local_kernel = numerics.local_boundary_kernels[type(bc)]
            # for edge in p.mesh.edges_on_region(region):
            for edge in mesh.edges_on(region):
                T, _ = edge["triangles"]
                for i in basis.dofs_on_element(T):
                    d_psi = basis.global_direction(i)
                    value = local_kernel.RHS(Edge(edge["M"], edge["l"], edge["N"], edge["T"]), d_psi, p.k)
                    rows.append(i)
                    values.append(value)

        b = np.zeros((basis.N_DOF,), dtype=np.complex128)
        np.add.at(b, rows, values)
        return b
