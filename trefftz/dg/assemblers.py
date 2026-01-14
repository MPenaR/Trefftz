'''Module for defininig the assemblers (if more than one)'''
from trefftz.mesh import TrefftzMesh
import numpy as np
# from trefftz.dg.serial_fluxes import SoundHard
from trefftz.mesh import EdgeType
from trefftz.dg.fluxes import FluxType, FluxKernels
from trefftz.dg.basis import TrefftzBasis

def serial_assembler(mesh: TrefftzMesh, N_theta: int):
    DOF = N_theta*mesh.n_triangles
    A = np.zeros((DOF, DOF), dtype=np.complex64)
    edges = mesh.edges
    # for edge in edges[edges["type"] == EdgeType.INNER]:
    #     for i in range(N_theta):
    #         for j in range(N_theta):
    #             A[i,j]

    # for edge in edges[edges["type"] == EdgeType.BOUNDARY]:
    #     if edge["flux_type"]==FluxType.SOUNDHARD:
    #         kernel = FluxKernels[edge["flux_type"]]
    #         f = kernel()
    for edge in edges:
        match edge["flux_type"]:

            case FluxType.SOUNDHARD:
                t_ID = edge["triangles"][0]


            case _:
                pass




# def sparse_assembler(mesh: Mesh):
#     for bc_types in mesh... 
        
