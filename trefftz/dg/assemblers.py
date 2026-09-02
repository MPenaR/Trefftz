from trefftz.mesh import TrefftzMesh
from scipy.sparse import coo_array
from trefftz.dg.serial_fluxes import SIGN as SSIGN
from trefftz.dg.block_fluxes import SIGN as BSIGN
import numpy as np
from trefftz.numpy_types import complex_array, int_array
from trefftz.dg.serial_fluxes import SerialLocalKernel, SerialNonLocalKernel, SerialTransmissionKernel
from trefftz.dg.block_fluxes import BlockLocalKernel, BlockNonLocalKernel, BlockTransmissionKernel
from trefftz.dg.basis import PlanewaveBasis
from typing import Protocol, Mapping, Any, Sequence
from dataclasses import dataclass
from trefftz.dg.boundary_conditions import BoundaryCondition
from abc import ABC, abstractmethod
from enum import IntEnum, StrEnum
from tqdm import tqdm


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

class Assembler[B: (StrEnum, IntEnum), num: Numerics](ABC):
    def __init__(self,
                 mesh: TrefftzMesh[B, Any], 
                 boundary_conditions: Mapping[B, BoundaryCondition],
                 numerics: num,
                 basis: PlanewaveBasis,
                 verbose: bool = True):
        
        self._mesh = mesh
        self._boundary_conditions = boundary_conditions
        self._numerics = numerics
        self._basis = basis
        self._regions_local_kernel = [region for region, bc in boundary_conditions.items() if type(bc) in numerics.local_boundary_kernels]
        self._regions_nonlocal_kernel = [region for region, bc in boundary_conditions.items() if type(bc) in numerics.nonlocal_boundary_kernels]
        self._regions_RHS_term = [region for region, bc in boundary_conditions.items() if bc.data is not None]
        self.verbose = verbose

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
        verbose = True
        for edge in tqdm(edges_on_region,
                         disable=not verbose,
                         unit="edge"):
            T, _ = edge["triangles"]
            for i in basis.dofs_on_element(T):
                for j in basis.dofs_on_element(T):
                    d_psi = basis.global_direction(i)
                    d_phi = basis.global_direction(j)
                    value = kernel.LHS(edge=edge, d_u=d_phi, d_v=d_psi, k=basis.k)
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

        verbose = self.verbose

        # boundary conditions implemented as local operators
        if verbose: 
            print('assembling local operators')
        for region in regions_local_kernel:
            if verbose: 
                print(f'region: {region}')
            edges_on_region = mesh.edges_on_boundary(region)
            bc = boundary_conditions[region]
            kernel = numerics.local_boundary_kernels[type(bc)]
            self.assemble_local_bc(edges_on_region, kernel, basis, rows, cols, values)

        refractive_index = basis.refractive_index

        # interior edges
        interior_kernel = numerics.interior_kernel

        for edge in tqdm(mesh.interior_edges,
                         desc="Interior",
                         disable=not verbose,
                         unit="edge"):
            for (i_v, T_v) in enumerate(edge["triangles"]):
                n_v = refractive_index.at(T_v)
                for (i_u, T_u) in enumerate(edge["triangles"]):
                    n_u = refractive_index.at(T_u)
                    sign = SSIGN((i_u, i_v))
                    for i in basis.dofs_on_element(T_v):
                        d_psi = basis.global_direction(i)
                        for j in basis.dofs_on_element(T_u):
                            d_phi = basis.global_direction(j)
                            val = interior_kernel.LHS(edge=edge, d_u=d_phi, d_v=d_psi, k=basis.k, sign=sign)
                            rows.append(i)
                            cols.append(j)
                            values.append(val)


        print("assembling non-local operators")
        # boundary conditions implemented as non-local operators
        for region in regions_nonlocal_kernel:
            bc = boundary_conditions[region]
            non_local_kernel = numerics.nonlocal_boundary_kernels[type(bc)]
            for edge_1 in tqdm(mesh.edges_on_boundary(region),
                               desc=f"NtD, {region}",
                               disable=not verbose,
                               unit="edge"):
                T_1, _ = edge_1["triangles"]
                for edge_2 in mesh.edges_on_boundary(region):
                    T_2, _ = edge_2["triangles"]
                    for j in basis.dofs_on_element(T_1):
                        for i in basis.dofs_on_element(T_2):
                            d_phi = basis.global_direction(j)
                            d_psi = basis.global_direction(i)
                            value = non_local_kernel.LHS(edge_u=edge_1,                                                                                edge_v=edge_2,                       # ... edge_u = edge_1, pero:
                                                         d_u=d_phi, d_v=d_psi, k=basis.k)     
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
            for edge in mesh.edges_on_boundary(region):
                T, _ = edge["triangles"]
                for i in basis.dofs_on_element(T):
                    d_psi = basis.global_direction(i)
                    value = local_kernel.RHS(edge=edge, d_v=d_psi, k=basis.k)
                    rows.append(i)
                    values.append(value)

        b = np.zeros((basis.N_DOF,), dtype=np.complex128)
        np.add.at(b, rows, values)
        return b

class BlockAssembler(Assembler[Any, BlockNumerics]):

    def assemble_LHS(self) -> coo_array:

        regions_local_kernel = self._regions_local_kernel
        regions_nonlocal_kernel = self._regions_nonlocal_kernel
        mesh = self._mesh
        boundary_conditions = self._boundary_conditions
        numerics = self._numerics
        basis = self._basis


        rows_dof: list[int_array] = []
        cols_dof: list[int_array] = []
        blocks: list[complex_array] = []
        refractive_index = basis.refractive_index

        verbose = self.verbose

        # boundary conditions implemented as local operators
        if verbose: 
            print('assembling local operators')
        for region in regions_local_kernel:
            if verbose: 
                print(f'region: {region}')
            edges_on_region = mesh.edges_on_boundary(region)
            bc = boundary_conditions[region]
            local_flux = numerics.local_boundary_kernels[type(bc)]

            for edge in tqdm(edges_on_region,
                            disable=not verbose,
                            unit="edge"):
                T, _ = edge["triangles"]
                dof = basis.dofs_on_element(T)
                D =  basis.global_direction(dof)
                n = refractive_index.at(T)
                block = local_flux.LHS(edge=edge, D=D, n=n, k=basis.k)
                rows_dof.append(dof)
                cols_dof.append(dof)
                blocks.append(block)


        # interior edges
        interior_flux = numerics.interior_kernel
        for edge in tqdm(mesh.interior_edges,
                         desc="Interior",
                         disable=not verbose,
                         unit="edge"):
            for (i_v, T_v) in enumerate(edge["triangles"]):
                for (i_u, T_u) in enumerate(edge["triangles"]):
                    sign = BSIGN((i_u, i_v))
                    u_dof = basis.dofs_on_element(T_u)
                    v_dof = basis.dofs_on_element(T_v)
                    D_v =  basis.global_direction(v_dof)
                    D_u =  basis.global_direction(u_dof)
                    n_u = refractive_index.at(T_u)
                    n_v = refractive_index.at(T_v)
                    block = interior_flux.LHS(edge=edge, D_u=D_u, n_u=n_u, D_v=D_v, n_v=n_v, k=basis.k, sign=sign)
                    rows_dof.append(v_dof)
                    cols_dof.append(u_dof)
                    blocks.append(block)


        print("assembling non-local operators")
        # boundary conditions implemented as non-local operators
        for region in regions_nonlocal_kernel:
            bc = boundary_conditions[region]
            non_local_flux = numerics.nonlocal_boundary_kernels[type(bc)]
            for edge_u in tqdm(mesh.edges_on_boundary(region),
                               desc=f"NtD, {region}",
                               disable=not verbose,
                               unit="edge"):
                T_u, _ = edge_u["triangles"]
                u_dof = basis.dofs_on_element(T_u)
                D_u =  basis.global_direction(u_dof)
                n_u = refractive_index.at(T_u)
                for edge_v in mesh.edges_on_boundary(region):
                    T_v, _ = edge_v["triangles"]
                    v_dof = basis.dofs_on_element(T_v)
                    D_v =  basis.global_direction(v_dof)
                    n_v = refractive_index.at(T_v)
                    block = non_local_flux.LHS(edge_u=edge_u, edge_v=edge_v, D_u=D_u, n_u=n_u, D_v=D_v, n_v=n_v, k=basis.k)
                    rows_dof.append(v_dof)
                    cols_dof.append(u_dof)
                    blocks.append(block)
        rows: list[int_array] = []
        cols: list[int_array] = []
        values: list[complex_array] = []

        for r, c, block in zip(rows_dof, cols_dof, blocks):
            rows.append(np.repeat(r, len(c)))  # 1 1 1 2 2 2 3 3 3...
            cols.append(np.tile(c, len(r)))  # 1 2 3 1 2 3 1 2 3...
            values.append(np.ravel(block, order="C"))  #C mayor order, i.e. column changes the fastest

        return coo_array((np.concatenate(values), (np.concatenate(rows), np.concatenate(cols))), shape=(basis.N_DOF, basis.N_DOF))

    def assemble_RHS(self) -> complex_array:

        rows_dof: list[int_array] = []
        blocks: list[complex_array] = []

        regions_RHS_term = self._regions_RHS_term
        mesh = self._mesh
        boundary_conditions = self._boundary_conditions
        numerics = self._numerics
        basis = self._basis

        for region in regions_RHS_term:  # I should check redefining this lists as sets or something like that, because of the local AND RHS
            bc = boundary_conditions[region]
            flux = numerics.local_boundary_kernels[type(bc)]
            for edge in mesh.edges_on_boundary(region):
                T, _ = edge["triangles"]
                v_dof = basis.dofs_on_element(T)
                rows_dof.append(v_dof)
                D_v = basis.global_direction(v_dof)
                block = flux.RHS(edge=edge, D=D_v, n=1., k=basis.k)
                blocks.append(block)
            
        rows = np.concatenate(rows_dof)
        values = np.concatenate(blocks)
        b = np.zeros((basis.N_DOF,), dtype=np.complex128)
        np.add.at(b, rows, values)
        return b

# class SerialFlux(ABC):
    
#     @abstractmethod
#     def LHS(self) -> complex:
#         pass

#     @abstractmethod
#     def RHS(self) -> complex:
#         pass


# class BlockFlux(ABC):
#     pass



# class GeneralAssembler[B]:

#     def __init__(self,
#                 mesh: TrefftzMesh[B, Any], 
#                 boundary_conditions: Mapping[B, BoundaryCondition],
#                 numerics: Numerics,
#                 basis: PlanewaveBasis,
#                 verbose: bool = True):
        
#         self._mesh = mesh
#         self._boundary_conditions = boundary_conditions
#         self._numerics = numerics
#         self._basis = basis
#         self._regions_local_kernel = [region for region, bc in boundary_conditions.items() if type(bc) in numerics.local_boundary_kernels]
#         self._regions_nonlocal_kernel = [region for region, bc in boundary_conditions.items() if type(bc) in numerics.nonlocal_boundary_kernels]
#         self._regions_RHS_term = [region for region, bc in boundary_conditions.items() if bc.data is not None]
#         self.verbose = verbose
#         self._rows: list[int] = []
#         self._columns: list[int] = []
#         self._values: list[complex] = []

#     # def assemble_block(self, flux: SerialFlux | BlockFlux, edges,
#     #                    dof_u: Sequence[int], dof_v: Sequence[int],
#     #                    D_u: float_array, D_v: float_array):
#     #     basis = self._basis
#     #     refractive_index = basis.refractive_index
#     #     match flux:
#     #         case SerialFlux():
#     #             for edge in tqdm(edges,
#     #                             desc=description,
#     #                             disable=not verbose,
#     #                             unit="edge"):
#     #                 for (i_v, T_v) in enumerate(edge["triangles"]):
#     #                     n_v = refractive_index.at(T_v)
#     #                     for (i_u, T_u) in enumerate(edge["triangles"]):
#     #                         n_u = refractive_index.at(T_u)
#     #                         sign = SSIGN((i_u, i_v))
#     #                         for i, d_v in zip(dof_v, D_v):
#     #                             for j, d_u in dof_u:
#     #                                 val = flux.LHS(edge=edge, d_u=d_u, n_u=n_u, d_v=d_v, n_v=n_v, k=basis.k, sign=sign)
#     #                                 self._rows.append(i)
#     #                                 self._columns.append(j)
#     #                                 self._values.append(val)
#     #         case BlockFlux():
#     #             for edge in tqdm(edges,
#     #                              desc=description,
#     #                 disable=not verbose,
#     #                 unit="edge"):
#     #                 for (i_v, T_v) in enumerate(edge["triangles"]):
#     #                     for (i_u, T_u) in enumerate(edge["triangles"]):
#     #                         sign = BSIGN((i_u, i_v))
#     #                         u_dof = basis.dofs_on_element(T_u)
#     #                         v_dof = basis.dofs_on_element(T_v)
#     #                         D_v =  basis.global_direction(v_dof)
#     #                         D_u =  basis.global_direction(u_dof)
#     #                         n_u = refractive_index.at(T_u)
#     #                         n_v = refractive_index.at(T_v)
#     #                         block = interior_flux.LHS(edge=edge, D_u=D_u, n_u=n_u, D_v=D_v, n_v=n_v, k=basis.k, sign=sign)
#     #                         rows_dof = np.repeat(v_dof, len(u_dof))
#     #                         cols_dof = np.tile(u_dof, len(v_dof))
#     #                         rows.extend(rows_dof)
#     #                         cols.extend(cols_dof)
#     #                         values.extend(block.ravel())
#     #         case _:
#     #             raise TypeError(
#     #                 f"Invalid flux type {type(interior_flux).__name__}: "
#     #                 "flux must inherit from SerialFlux or BlockFlux.")                



#     def assemble_LHS(self) -> coo_array:

#         boundaries_with_local_flux = self._regions_local_kernel
#         boundaries_with_nonlocal_flux = self._regions_nonlocal_kernel
#         mesh = self._mesh
#         boundary_conditions = self._boundary_conditions
#         numerics = self._numerics
#         basis = self._basis
#         refractive_index = basis.refractive_index


#         rows: list[int] = []
#         cols: list[int] = []
#         values: list[complex] = []

#         verbose = self.verbose

#         # boundary conditions implemented as local operators
#         if verbose: 
#             print('assembling local operators')

#         for boundary in boundaries_with_local_flux:            
#             edges_on_boundary = mesh.edges_on_boundary(boundary)
#             bc = boundary_conditions[boundary]
#             local_flux = numerics.local_boundary_kernels[type(bc)]
#             if verbose: 
#                 print(f'Boundary: {boundary},\n Boundary condition: {bc},\n Flux implementation {local_flux}')

#             match local_flux: 
#                 case SerialFlux():
#                     for edge in tqdm(edges_on_boundary,
#                          disable=not verbose,
#                          unit="edge"):
#                         T, _ = edge["triangles"]
#                         dof = basis.dofs_on_element(T)
#                         for i in dof:
#                             for j in dof:
#                                 d_psi = basis.global_direction(i)
#                                 d_phi = basis.global_direction(j)
#                                 value = local_flux.LHS(edge=edge, d_u=d_phi, d_v=d_psi, k=basis.k)
#                                 rows.append(i)
#                                 cols.append(j)
#                                 values.append(value)
                
#                 case BlockFlux():
#                     for edge in tqdm(edges_on_boundary,
#                             disable=not verbose,
#                             unit="edge"):
#                         T, _ = edge["triangles"]
#                         dof = basis.dofs_on_element(T)
#                         D =  basis.global_direction(dof)
#                         n = refractive_index.at(T)
#                         block = local_flux.LHS(edge=edge, D=D, n=n, k=basis.k)
#                         rows_dof = np.repeat(dof, len(dof))
#                         cols_dof = np.tile(dof, len(dof))
#                         rows.extend(rows_dof)
#                         cols.extend(cols_dof)
#                         values.extend(block.ravel())
#                 case _:
#                     raise TypeError(
#                         f"Invalid flux type {type(local_flux).__name__}: "
#                         "flux must inherit from SerialFlux or BlockFlux.")                
              

#         # interior edges
#         interior_edges = mesh.interior_edges        
#         interior_flux = numerics.interior_kernel
#         if verbose: 
#             print(f'Interior edges,\n Flux implementation {interior_flux}')
        
#         match interior_flux:
#             case SerialFlux():
#                 for edge in tqdm(interior_edges,
#                                 desc="Interior",
#                                 disable=not verbose,
#                                 unit="edge"):
#                     for (i_v, T_v) in enumerate(edge["triangles"]):
#                         n_v = refractive_index.at(T_v)
#                         for (i_u, T_u) in enumerate(edge["triangles"]):
#                             n_u = refractive_index.at(T_u)
#                             sign = SSIGN((i_u, i_v))
#                             for i in basis.dofs_on_element(T_v):
#                                 d_psi = basis.global_direction(i)
#                                 for j in basis.dofs_on_element(T_u):
#                                     d_phi = basis.global_direction(j)
#                                     val = interior_flux.LHS(edge=edge, d_u=d_phi, n_u=n_u, d_v=d_psi, n_v=n_v, k=basis.k, sign=sign)
#                                     rows.append(i)
#                                     cols.append(j)
#                                     values.append(val)
#             case BlockFlux():
#                 for edge in tqdm(mesh.interior_edges,
#                     desc="Interior",
#                     disable=not verbose,
#                     unit="edge"):
#                     for (i_v, T_v) in enumerate(edge["triangles"]):
#                         for (i_u, T_u) in enumerate(edge["triangles"]):
#                             sign = BSIGN((i_u, i_v))
#                             u_dof = basis.dofs_on_element(T_u)
#                             v_dof = basis.dofs_on_element(T_v)
#                             D_v =  basis.global_direction(v_dof)
#                             D_u =  basis.global_direction(u_dof)
#                             n_u = refractive_index.at(T_u)
#                             n_v = refractive_index.at(T_v)
#                             block = interior_flux.LHS(edge=edge, D_u=D_u, n_u=n_u, D_v=D_v, n_v=n_v, k=basis.k, sign=sign)
#                             rows_dof = np.repeat(v_dof, len(u_dof))
#                             cols_dof = np.tile(u_dof, len(v_dof))
#                             rows.extend(rows_dof)
#                             cols.extend(cols_dof)
#                             values.extend(block.ravel())
#             case _:
#                 raise TypeError(
#                     f"Invalid flux type {type(interior_flux).__name__}: "
#                     "flux must inherit from SerialFlux or BlockFlux.")                

#         print("assembling non-local operators")
#         # boundary conditions implemented as non-local operators
#         for boundary in boundaries_with_nonlocal_flux:
#             bc = boundary_conditions[boundary]
#             non_local_flux = numerics.nonlocal_boundary_kernels[type(bc)]
#             if verbose: 
#                 print(f'Boundary: {boundary},\n Boundary condition: {bc},\n Flux implementation {non_local_flux}')
            
#             match non_local_flux:
#                 case SerialFlux():
#                     for edge_1 in tqdm(mesh.edges_on_boundary(boundary),
#                                     desc=f"NtD, {boundary}",
#                                     disable=not verbose,
#                                     unit="edge"):
#                         T_1, _ = edge_1["triangles"]
#                         for edge_2 in mesh.edges_on_boundary(boundary):
#                             T_2, _ = edge_2["triangles"]
#                             for j in basis.dofs_on_element(T_1):
#                                 for i in basis.dofs_on_element(T_2):
#                                     d_phi = basis.global_direction(j)
#                                     d_psi = basis.global_direction(i)
#                                     value = non_local_flux.LHS(edge_u=edge_1,                                                                                edge_v=edge_2,                       # ... edge_u = edge_1, pero:
#                                                                 d_u=d_phi, d_v=d_psi, k=basis.k)     
#                                     rows.append(i)                                                    
#                                     cols.append(j)                                                    
#                                     values.append(value)
#                 case BlockFlux():
#                     for edge_u in tqdm(mesh.edges_on_boundary(boundary),
#                                         desc=f"NtD, {boundary}",
#                                         disable=not verbose,
#                                         unit="edge"):
#                         T_u, _ = edge_u["triangles"]
#                         u_dof = basis.dofs_on_element(T_u)
#                         D_u =  basis.global_direction(u_dof)
#                         n_u = refractive_index.at(T_u)
#                         for edge_v in mesh.edges_on_boundary(boundary):
#                             T_v, _ = edge_v["triangles"]
#                             v_dof = basis.dofs_on_element(T_v)
#                             D_v =  basis.global_direction(v_dof)
#                             n_v = refractive_index.at(T_v)
#                             block = non_local_flux.LHS(edge_u=edge_u, edge_v=edge_v, D_u=D_u, n_u=n_u, D_v=D_v, n_v=n_v, k=basis.k)
#                             rows_dof = np.repeat(v_dof, len(u_dof))
#                             cols_dof = np.tile(u_dof, len(v_dof))
#                             rows.extend(rows_dof)
#                             cols.extend(cols_dof)
#                             values.extend(block.ravel())
#                 case _:
#                     raise TypeError(
#                     f"Invalid flux type {type(non_local_flux).__name__}: "
#                     "flux must inherit from SerialFlux or BlockFlux.")                



#         return coo_array((np.asarray(values), (np.asarray(rows), np.asarray(cols))), shape=(basis.N_DOF, basis.N_DOF))

#     def assemble_RHS(self) -> complex_array:

#         rows: list[int] = []
#         values: list[complex] = []

#         regions_RHS_term = self._regions_RHS_term
#         mesh = self._mesh
#         boundary_conditions = self._boundary_conditions
#         numerics = self._numerics
#         basis = self._basis

#         for region in regions_RHS_term:  # I should check redefining this lists as sets or something like that, because of the local AND RHS
#             bc = boundary_conditions[region]
#             local_kernel = numerics.local_boundary_kernels[type(bc)]
#             for edge in mesh.edges_on_boundary(region):
#                 T, _ = edge["triangles"]
#                 for i in basis.dofs_on_element(T):
#                     d_psi = basis.global_direction(i)
#                     value = local_kernel.RHS(edge=edge, d_v=d_psi, k=basis.k)
#                     rows.append(i)
#                     values.append(value)

#         b = np.zeros((basis.N_DOF,), dtype=np.complex128)
#         np.add.at(b, rows, values)
#         return b
