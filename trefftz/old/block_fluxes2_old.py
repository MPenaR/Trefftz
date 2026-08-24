r"""
module for implementing blocks of fluxes, i.e. matrices containing all the fluxes of the same type for a 
given pair of cells.
"""

from trefftz.numpy_types import float_array, complex_array
from numpy import dot, exp, sqrt, sinc, subtract, add, pi, outer
import numpy as np
from dataclasses import dataclass
from trefftz.dg.basis import PlanewaveBasis

@dataclass(slots=True)
class Edge:
    P: float_array
    Q: float_array
    l: float
    T: float
    N: float
    M: float


def from_edge_to_Edge( edge: np.void) -> Edge:
    return Edge(edge["P"], edge["Q"], edge["l"], edge["T"], edge["N"], edge["M"])



def SoundHard_block_planewave(k: complex, l: float, N: float_array, T: float_array, M: float_array, d: float_array, d_d: float_array, d_1: float) -> complex_array:
    
    r"""
    Computes the block for an edge in a sound hard boundary.
    That is it computes the matrix :math:`\mathbf{M}=(M_{mn})` with:    
    
    .. math::
    
        M_{mn}=\boxed{(-ikl\left(1+d_{1}\mathbf{d}_{n}\cdot\mathbf{n}\right)\mathbf{d}_{m}\cdot\mathbf{n}e^{ik\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\cdot\mathbf{M}}\mathrm{sinc}\left(\frac{kl}{2\pi}\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\cdot\boldsymbol{\tau}\right)}


    Parameters:
    -----------

    - k : complex
        Wavenumber
    - l: float
        Length of the edge
    - N: float_array
        Normal of the edge
    - T: float_array
        Tangent of the edge
    - M: float_array
        Midpoint of the edge
    - d : float_array
        Set of directions
    - d_d : float_array
        Nd x Nd x 2 "Matrix" of differences of directions.
    - d_1 : float
        Stabilyzing parameter.
    """

    return -1j*k*l*outer(dot(d, N), (1 + d_1*dot(d, N)))*exp(1j*k*dot(d_d, M))*sinc(k*l/(2*pi)*dot(d_d, T))

def SoundHard_block(edge: Edge, basis: PlanewaveBasis, d_1: float) -> complex_array:
    
    r"""
    Computes the block for an edge in a sound hard boundary assuming a plane wave basis.
    That is it computes the matrix :math:`\mathbf{M}=(M_{mn})` with:    
    
    .. math::
    
        M_{mn}=\boxed{-ikl\left(1+d_{1}\mathbf{d}_{n}\cdot\mathbf{n}\right)\mathbf{d}_{m}\cdot\mathbf{n}e^{ik\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\cdot\mathbf{M}}\mathrm{sinc}\left(\frac{kl}{2\pi}\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\cdot\boldsymbol{\tau}\right)}


    Parameters:
    -----------

    - k : complex
        Wavenumber
    - edge : Edge
        Edge
    - d : float_array
        Set of directions
    - d_d : float_array
        Nd x Nd x 2 "Matrix" of differences of directions.
    - d_1 : float
        Stabilyzing parameter.
    """


    return SoundHard_block_planewave(k=basis.k,  l=edge.l, N=edge.N, T=edge.T, M=edge.M, d=basis.D, d_d=basis.D_D, d_1=d_1)
    

# def Inner_block(k : complex, edge : Edge, d : float_array, a : float, b : float, n_A : float, n_B : float) -> complex_array:
#     r"""
#     Computes the block for an inner edge.

#     Parameters:
#     -----------

#     - k : complex
#         Wavenumber
#     - edge : Edge
#         Edge
#     - d : float_array
#         Set of directions
#     - a : float
#         Stabilyzing parameter.    
#     - b : float
#         Stabilyzing parameter.    
# """
#     l = edge.l
#     N = edge.N
#     M = edge.M
#     T = edge.T

#     #change later
#     n_m = 1
#     n_n = 1
    
#     I = -1j*k*l*(a + add.outer(sqrt(n_m)*dot(d,N),sqrt(n_n)*dot(d,N))/2 + b*outer(sqrt(n_m)*dot(d,N),sqrt(n_n)*dot(d,N))) \
#     *exp(-1j*k* subtract.outer(sqrt(n_m)*dot(d,M),sqrt(n_n)*dot(d,M)))                                             \
#     *sinc(l*k/(2*pi)*subtract.outer(sqrt(n_m)*dot(d,T),sqrt(n_n)*dot(d,T)))

#     return I

# def Radiating_local_block(k : complex, edge : Edge, d : float_array, d_d : float_array, d_2 : float) -> complex_array:
#     r"""
#     Computes the same triangle block for a radiating edge.

#     Parameters:
#     -----------

#     - k : complex
#         Wavenumber
#     - edge : Edge
#         Edge
#     - d : float_array
#         Set of directions
#     - d_d : float_array
#         Nd x Nd x 2 "Matrix" of differences of directions.
#     - d_2 : float
#         Stabilyzing parameter.    
# """

#     l = edge.l
#     N = edge.N
#     M = edge.M
#     T = edge.T



#     I = -1j*k*l*(d_2 + dot(d, N))*exp(1j*k*dot(d_d, M))*sinc(k*l/(2*pi)*dot(d_d, T))
#     return I
