'''Module for defininig the assemblers (if more than one)'''
from trefftz.mesh import Mesh
import numpy as np
from trefftz.fem.serial_fluxes import Inner



def serial_assembler(mesh: Mesh, N_theta: int):
    DOF = N_theta*mesh.n_triangles
    A = np.zeros((DOF, DOF), dtype=np.complex64)
    edges = mesh.edges
    for edge in edges[np.logical_not(edges["boundary"])]:






# def sparse_assembler(mesh: Mesh):
#     for bc_types in mesh... 
        
