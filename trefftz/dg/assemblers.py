'''Module for defininig the assemblers (if more than one)'''
from trefftz.mesh import TrefftzMesh
from trefftz.numpy_types import complex_array, float_array
import numpy as np
# from trefftz.dg.serial_fluxes import SoundHard
from trefftz.mesh import EdgeType
from trefftz.dg.fluxes import FluxType, FluxKernels, SoundHard
from trefftz.dg.basis import TrefftzBasis
from typing import NamedTuple


class BasisFunction(NamedTuple):
    d: float_array
    n: float


# BasisFunction = namedtuple("BasisFunction", ["d", "n"])

def serial_assembler(mesh: TrefftzMesh, basis: TrefftzBasis) -> complex_array:
    N_DOF = basis.N_DOF
    A = np.zeros((N_DOF, N_DOF), dtype=np.complex64)
    edges = mesh.edges
    k = basis.k
    # for edge in edges[edges["type"] == EdgeType.INNER]:
    #     for i in range(N_theta):
    #         for j in range(N_theta):
    #             A[i,j]

    # for edge in edges[edges["type"] == EdgeType.BOUNDARY]:
    #     if edge["flux_type"]==FluxType.SOUNDHARD:
    #         kernel = FluxKernels[edge["flux_type"]]
    #         f = kernel()
    D = basis.D
    for edge in edges:
        match edge["flux_type"]:
            case FluxType.SOUNDHARD:
                T = edge["triangles"][0]
                for i in basis.dofs_on_triangle(T):
                    for j in basis.dofs_on_triangle(T):
                        psi = BasisFunction(d=D[j, :], n=1)
                        phi = BasisFunction(d=D[i, :], n=1)
                        A[i, j] = SoundHard(phi, psi, k, edge, d_1=0.5)
            # case FluxType.TRANSMISSION:
            #     T_plus, T_minus = edge["triangles"]
            #     for i in basis.dofs_on_triangle(T):
            #         for j in basis.dofs_on_triangle(T):
            #             psi = BasisFunction(d=D[j, :], n=1)
            #             phi = BasisFunction(d=D[i, :], n=1)
            #             A[i, j] = SoundHard(phi, psi, k, edge, d_1=0.5)
            case _:
                pass
    return A




# def sparse_assembler(mesh: Mesh):
#     for bc_types in mesh... 
        
