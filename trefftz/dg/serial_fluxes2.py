r"""
Module for the evaluation of fluxes for a single pair of test and trial functions. It was the initial
implementation of the code and now serves as a test for the vectorized implementation. It is tested
against numerical computations of the fluxes.

All around this page, :math:`\Phi_{n}(\mathbf{x})=\exp(ik\mathbf{x}\cdot\mathbf{d}_n)` and :math:`\Psi_{m}(\mathbf{x})=\exp(ik\mathbf{x}\cdot\mathbf{d}_m)` 
are the trial and test functions. :math:`l` is the length of an edge, :math:`\mathbf{M}` its midpoint and :math:`\boldsymbol{\tau}` and :math:`\mathbf{n}` 
its tangent and normal unitary vectors. 
"""


import numpy as np
from numpy import dot, sinc, pi, exp, sqrt, conj
from dataclasses import dataclass
from trefftz.numpy_types import float_array
from typing import Mapping


# @dataclass
# class Function:
#     d: np.ndarray[tuple[int], np.dtype[np.floating]]
#     n: float


def SoundHard(d_m, d_n, k: float, edge, flux_parameters: Mapping[str, float] = {"d_1": 0.5}) -> complex:
    r"""
    Computes the flux on a sound-hard boundary, that is:

    .. math::
    
        \int_{E}\left(\varphi_n(\mathbf{x})+\frac{d_{1}}{ik}\nabla \varphi_n(\mathbf{x})\cdot\mathbf{n}\right)\overline{\nabla \psi_m(\mathbf{x})\cdot\mathbf{n}}\,\mathrm{d}S_{\mathbf{x}}

    This quantity can be exactly evaluated as:

    .. math::
    
        \boxed{-ikl\left(1+d_{1}\mathbf{d}_{n}\cdot\mathbf{n}\right)\mathbf{d}_{m}\cdot\mathbf{n}e^{ik\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\cdot\mathbf{M}}\mathrm{sinc}\left(\frac{kl}{2\pi}\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\cdot\boldsymbol{\tau}\right)}


    Parameters
    ----------
    phi : Function
        Trial function.
    psi : Function
        Test function.
    k : float
        Wavenumber.
    edge : Edge
        Edge parameters.
    d_1 : float
        Stabilyzing parameter.

    Returns
    -------
    I : complex
        The integral.
    
    """
    
    M = edge["M"]
    l = edge["l"]
    N = edge["N"]
    T = edge["T"]

    d_1 = flux_parameters["d_1"]

    return -1j*k*l*(1 + d_1 * dot(d_n, N))*dot(d_m, N)*exp(1j*k*dot(d_n - d_m, M)) * sinc(k*l/(2*pi)*dot(d_n-d_m, T))
    



def InnerEasy(d_m, d_n, k: float, edge, flux_parameters: Mapping[str, float] = {"a": 0.5, "b": 0.5}) -> complex:
    r"""
    Computes the flux on a inner facet with respect to the degrees
    of freedom from the same cell, that is:
    
    .. math::
        \int_E \left(\left(\varphi_n(\mathbf{x})+\frac{b}{ik}\nabla\varphi_n(\mathbf{x})\cdot\mathbf{n}\right)\overline{\nabla\psi_m(\mathbf{x})\cdot\mathbf{n}}- \left(\vphantom{\frac{1}{2}}aik\varphi_n(\mathbf{x})+\nabla\varphi_n(\mathbf{x})\cdot\mathbf{n}\right)\overline{\psi_m(\mathbf{x})}\right) \,\mathrm{d}S_\mathbf{x}    

    Parameters
    ----------
    phi : Function
        Trial function.
    psi : Function
        Test function.
    k : float
        Wavenumber.
    edge : Edge
        Edge parameters.
    a : float
        Stabilyzing parameter.
    b : float
        Stabilyzing parameter.

    Returns
    -------
    I : complex
        The integral.


    """

    k_m = k
    k_n = k

    a = flux_parameters["a"]
    b = flux_parameters["b"]

    
    M = edge["M"]
    l = edge["l"]
    N = edge["N"]
    T = edge["T"]

    # I = -1j*l/2*(2*a*k + k_n*dot(d_n, N) + k_m*dot(d_m, N) + 2*b/k*k_n*dot(d_n, N)*k_m*dot(d_m, N))*exp(1j*dot(k_n*d_n - k_m*d_m, M))*sinc(l/(2*pi)*dot(k_n*d_n - k_m*d_m,T))
    return -1j*k*l*( 1/2*( k_n/k*dot(d_n, N) + k_m/k*dot(d_m, N)) + a + b*k_n/k*k_m/k*dot(d_n, N)*dot(d_m, N))*exp(1j*dot(k_n*d_n - k_m*d_m, M))*sinc(l/(2*pi)*dot(k_n*d_n - k_m*d_m,T))

def Inner2(d_m, d_n, k: float, edge, flux_parameters: Mapping[str, float] = {"a": 0.5, "b": 0.5}) -> complex:
    r"""
    Computes the flux on a inner facet with respect to the degrees
    of freedom from the same cell, that is:
    
    .. math::
        \int_E \left(\left(\varphi_n(\mathbf{x})+\frac{b}{ik}\nabla\varphi_n(\mathbf{x})\cdot\mathbf{n}\right)\overline{\nabla\psi_m(\mathbf{x})\cdot\mathbf{n}}- \left(\vphantom{\frac{1}{2}}aik\varphi_n(\mathbf{x})+\nabla\varphi_n(\mathbf{x})\cdot\mathbf{n}\right)\overline{\psi_m(\mathbf{x})}\right) \,\mathrm{d}S_\mathbf{x}    

    Parameters
    ----------
    phi : Function
        Trial function.
    psi : Function
        Test function.
    k : float
        Wavenumber.
    edge : Edge
        Edge parameters.
    a : float
        Stabilyzing parameter.
    b : float
        Stabilyzing parameter.

    Returns
    -------
    I : complex
        The integral.


    """

    a = flux_parameters["a"]
    b = flux_parameters["b"]

    
    M = edge["M"]
    l = edge["l"]
    N = edge["N"]
    T = edge["T"]

    return -1j*k*l*((1/2+b*dot(d_n, N))*dot(d_m, N) + 1/2*dot(d_n, N) + a)*exp(1j*k*dot(d_n-d_m, M))*sinc(k*l/(2*pi)*dot(d_n-d_m, T))


# def Inner(d_m, d_n, k: float, edge, flux_parameters: Mapping[str, float] = {"a": 0.5, "b": 0.5}) -> complex:
#     r"""
#     Computes the flux on a inner facet with respect to the degrees
#     of freedom from the same cell, that is:
    
#     .. math::
#         \int_E \left(\left(\varphi_n(\mathbf{x})+\frac{b}{ik}\nabla\varphi_n(\mathbf{x})\cdot\mathbf{n}\right)\overline{\nabla\psi_m(\mathbf{x})\cdot\mathbf{n}}- \left(\vphantom{\frac{1}{2}}aik\varphi_n(\mathbf{x})+\nabla\varphi_n(\mathbf{x})\cdot\mathbf{n}\right)\overline{\psi_m(\mathbf{x})}\right) \,\mathrm{d}S_\mathbf{x}    

#     Parameters
#     ----------
#     phi : Function
#         Trial function.
#     psi : Function
#         Test function.
#     k : float
#         Wavenumber.
#     edge : Edge
#         Edge parameters.
#     a : float
#         Stabilyzing parameter.
#     b : float
#         Stabilyzing parameter.

#     Returns
#     -------
#     I : complex
#         The integral.


#     """

#     k_n = k * sqrt(phi.n)
#     k_m = k * sqrt(psi.n)


#     a = flux_parameters["a"]
#     b = flux_parameters["b"]

    
#     M = edge["M"]
#     l = edge["l"]
#     N = edge["N"]
#     T = edge["T"]

#     # I = -1j*l/2*(2*a*k + k_n*dot(d_n, N) + k_m*dot(d_m, N) + 2*b/k*k_n*dot(d_n, N)*k_m*dot(d_m, N))*exp(1j*dot(k_n*d_n - k_m*d_m, M))*sinc(l/(2*pi)*dot(k_n*d_n - k_m*d_m,T))
#     return -1j*k*l*( 1/2*( k_n/k*dot(d_n, N) + k_m/k*dot(d_m, N)) + a + b*k_n/k*k_m/k*dot(d_n, N)*dot(d_m, N))*exp(1j*dot(k_n*d_n - k_m*d_m, M))*sinc(l/(2*pi)*dot(k_n*d_n - k_m*d_m,T))


def Radiating_local(d_m, d_n, k: float, edge, flux_parameters: Mapping[str, float] = {"d_2": 0.5}) -> complex:
    r"""
    Computes the flux on a radiating boundary with respect to the degrees
    of freedom from the same cell, that is:

    TODO: it is assuming that the radiating boundary consists of a vertical segment. This should be easy to generalize.
    
    .. math::
    
        -\int_{E}\left(d_{2}ik\Phi_n(\mathbf{x})+\nabla \Phi_n(\mathbf{x})\cdot\mathbf{n}\right)\overline{\Psi_n(\mathbf{x})}\,\mathrm{d}S_{\mathbf{x}}

    which can be computed as:

    .. math::
    
        \boxed{-ikl\left(d_{2}+\mathbf{d}_{n}\cdot\mathbf{n}\right)e^{ik\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\cdot\mathbf{M}}\mathrm{sinc}\left(\frac{kl}{2\pi}\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\mathbf{j}\right)}

    Parameters
    ----------
    phi : Function
        Trial function.
    psi : Function
        Test function.
    k : float
        Wave number.
    edge : Edge
        Edge parameters.
    d_2 : float
        Stabilyzing parameter.

    Returns
    -------
    I : complex
        The integral.
    
    """
    d_2 = flux_parameters["d_2"]

    M = edge["M"]
    l = edge["l"]
    N = edge["N"]
    T = edge["T"]

    return -1j*k*l*(d_2 + dot(d_n, N))*exp(1j*k*dot(d_n - d_m, M))*sinc(k*l/(2*pi)*dot(d_n-d_m, T))


def Radiating_nonlocal(d_n, d_m, k : float, edge_u, edge_v,
                       NtD_modes: int, H: float = 1.0, flux_parameters: Mapping[str, float] = {"d_2": 0.5}) -> complex:
    r"""
    Computes the flux on a radiating boundary with respect to the degrees
    of freedom from another cell, that is:

    TODO: it is assuming that the radiating boundary consists of a vertical segment. This should be easy to generalize.
    
    Parameters
    ----------
    phi : Function
        Trial function.
    psi : Function
        Test function.
    k : float
        Wave number.
    edge_u : Edge
        Edge of the triangle associated to the trial function.
    edge_v : Edge
        Edge of the triangle associated to the test function.
    d_2 : float
        Stabilyzing parameter.
    NtD_modes : int
        Number of modes for the approximation of the NtD map.
    H : float
        height of the waveguide. 

    Returns
    -------
    I : complex
        The integral.

    
    """

    d_2 = flux_parameters["d_2"]

    M_u = edge_u["M"]
    l_u = edge_u["l"]
    # N = edge["N"]
    # T = edge["T"]

    M_v = edge_v["M"]
    l_v = edge_v["l"]
    # N = edge["N"]
    # T = edge["T"]


    N = edge_u["N"]
    T = edge_u["T"]

    I1 = -1j*k*H*dot(d_n,N)*dot(d_m,N)*d_2*exp(1j*k*(dot(d_n,M_u) - dot(d_m,M_v)))*l_u/H*l_v/H*(
        sinc(k*l_u/(2*pi)*d_n[1])*sinc(k*l_v/(2*pi)*d_m[1]) + 1/2*sum( [ k**2 / abs(sqrt(complex(k**2 - (s*pi/H)**2)))**2 * (
        exp( 1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] + s*l_u/(2*H)) + exp(-1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] - s*l_u/(2*H)) ) *(
        exp(-1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] + s*l_v/(2*H)) + exp( 1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] - s*l_v/(2*H)) )
        for s in range(1,NtD_modes)]) )
    
    I2 = -1j*k*H*dot(d_n,N)*(dot(d_m,N)-d_2)*exp(1j*k*(dot(d_n,M_u) - dot(d_m,M_v)))*l_u/H*l_v/H*(
        sinc(k*l_u/(2*pi)*d_n[1])*sinc(k*l_v/(2*pi)*d_m[1]) + 1/2*sum( [ k / sqrt(complex(k**2 - (s*pi/H)**2)) * (
        exp( 1j*s*pi*M_u[1]/H)*sinc(k*l_u/(2*pi)*d_n[1] + s*l_u/(2*H)) + exp(-1j*s*pi*M_u[1]/H)*sinc(k*l_u/(2*pi)*d_n[1] - s*l_u/(2*H)) ) *(
        exp(-1j*s*pi*M_v[1]/H)*sinc(k*l_v/(2*pi)*d_m[1] + s*l_v/(2*H)) + exp( 1j*s*pi*M_v[1]/H)*sinc(k*l_v/(2*pi)*d_m[1] - s*l_v/(2*H)) )
        for s in range(1,NtD_modes)]) )
    
    I3 = 1j*k*H*dot(d_m,N)*d_2*exp(1j*k*(dot(d_n,M_u) - dot(d_m,M_v)))*l_u/H*l_v/H*(
        sinc(k*l_u/(2*pi)*d_n[1])*sinc(k*l_v/(2*pi)*d_m[1]) + 1/2*sum( [ k / conj(sqrt(complex(k**2 - (s*pi/H)**2))) * (
        exp( 1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] + s*l_u/(2*H)) + exp(-1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] - s*l_u/(2*H)) ) *(
        exp(-1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] + s*l_v/(2*H)) + exp( 1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] - s*l_v/(2*H)) )
        for s in range(1,NtD_modes)]) )

    return  I1 + I2 + I3

def mode_RHS(d_m, edge, k, H, d_2, t):
    d = d_m
    d_y = d[1]
    N = edge["N"]
    M = edge["M"]
    l = edge["l"]


    beta = np.emath.sqrt(k**2 - (t*pi/H)**2)

    F = 1j*k*l*exp(1j*beta*M[0])*exp(-1j*k*dot(d, M))*(dot(d, N) - d_2)*(exp( 1j*pi*t/H*M[1])*sinc(t*l/(2*H) - k*l*d_y/(2*pi)) + 
                                                                         exp(-1j*pi*t/H*M[1])*sinc(t*l/(2*H) + k*l*d_y/(2*pi)))

    if t == 0:
        S = 2j*k*l*dot(d, N)*d_2*exp(1j*(beta*M[0]-k*dot(d, M)))*sinc(k*l/(2*pi)*d[1])

    else:
        S = 1j*k*l*dot(d, N)*d_2*exp(1j*(beta*M[0]-k*dot(d, M))) * (k/conj(beta) *
                                                                        (exp( 1j*t*pi*M[1]/H)*sinc(t*l/(2*H) - k*l*d[1]/(2*pi)) + 
                                                                         exp(-1j*t*pi*M[1]/H)*sinc(t*l/(2*H) + k*l*d[1]/(2*pi))))
 
    return F + S


def mode_RHS2(d_m, edge, k, H, d_2, t):
    d = d_m
    d_x, d_y = d
    N = edge["N"]
    M = edge["M"]
    M_x, M_y = M
    T = edge["T"]
    l = edge["l"]

    beta = np.emath.sqrt(k**2 - (t*pi/H)**2)

    if t == 0:
        I = 1/sqrt(H)*2*1j*k*l*(dot(d, N) - d_2*(1-dot(d, N)))*exp(1j*k*(M_x - dot(d, M)))*sinc(k*l/(2*pi)*dot(d, T))
    else:
        I = sqrt(2/H)*1j*k*l*(dot(d, N) - d_2*(1-k/conj(beta)*dot(d, N)))*exp(1j*(beta*M_x - k*dot(d, M)))*( exp( 1j*t*pi*dot(M,T)/H)*sinc(k*l/(2*pi)*dot(d, T) - t*(l/(2*H)))+
                                                                                                             exp(-1j*t*pi*dot(M,T)/H)*sinc(k*l/(2*pi)*dot(d, T) + t*(l/(2*H))))
    return I
    #     LOCAL = 1/sqrt(H)*2*1j*k*l*(dot(d, N) - d_2*(1-dot(d, N)))*exp(1j*k*(M_x - dot(d, M)))*sinc(k*l/(2*pi)*dot(d, T))
    #     NON_LOCAL = 1/sqrt(H)*2*1j*k*l**exp(1j*k*(M_x - dot(d, M)))*sinc(k*l/(2*pi)*dot(d, T))
    # else: 
    #     LOCAL = sqrt(2/H)*1j*k*l*(dot(d, N) - d_2)*exp(1j*(beta*M_x - k*dot(d, M)))*( exp( 1j*t*pi*M_y/H)*sinc(k*l/(2*pi)*dot(d, T) - t*(l/(2*H)))+
    #                                                                                   exp(-1j*t*pi*M_y/H)*sinc(k*l/(2*pi)*dot(d, T) + t*(l/(2*H))))
    #     NON_LOCAL = sqrt(2/H)*1j*k*l*d_2*dot(d, N)*exp(1j*(beta*M_x - k*dot(d, M)))*k/conj(beta)*( exp( 1j*t*pi*M_y/H)*sinc(k*l/(2*pi)*dot(d, T) - t*(l/(2*H)))+
    # return LOCAL + NON_LOCAL
