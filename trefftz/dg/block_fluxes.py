r"""
Module for the evaluation of fluxes in a vectorized manner over all the pair of test and trial functions corresponding to a given edge.
"""

from trefftz.numpy_types import float_array, complex_array
from enum import Enum
from typing import Protocol
from trefftz.dg.kernels.block import I_uv, I_duv, I_udv, I_dudv, I_pw_dv, I_pw_v


class SIGN(Enum):
    '''Sign for the transmission kernel
    where PP (plus plus) stands for both trial
    and test function coming from the plus triangle
    whereas PM stands for the test function (psi) coming
    from the + triangle and the trial function (phi) from
    the - one.'''
    PP = (0, 0)
    PM = (0, 1)
    MP = (1, 0)
    MM = (1, 1)


class BlockTransmissionKernel(Protocol):
    def LHS(self, edge, D_u: float_array, D_v: float_array, k: float, sign: SIGN) -> complex_array:
        ...


class BlockLocalKernel(Protocol):
    def LHS(self, edge, D_u: float_array, D_v: float_array, k: float) -> complex_array:
        ...

    def RHS(self, edge, D_v: float_array, k: float) -> complex_array:
        ...


class BlockNonLocalKernel(Protocol):
    def LHS(self, edge_u, edge_v, D_u: float_array, D_v: float_array, k: float) -> complex_array:
        ...


class NeumannFlux:
    r"""
    Computes the flux on a Neumann boundary condition, that is:

    .. math::
    
        \int_{E}\left(\varphi_n(\mathbf{x})+\frac{d_{1}}{ik}\nabla \varphi_n(\mathbf{x})\cdot\mathbf{n}\right)\overline{\nabla \psi_m(\mathbf{x})\cdot\mathbf{n}}\,\mathrm{d}S_{\mathbf{x}}

    Parameters
    ----------
    edge : Edge or Arc
        Edge parameters.
    D_u : (float, float) array
        Propatagion direction of the trial function.
    D_v : (float, float) array
        Propatagion direction of the test function.
    k : float
        Wavenumber.
    d_1 : float
        Stabilyzing parameter.

    Returns
    -------
    I : complex
        The integral.
    
    """
    def __init__(self, d_1: float):
        self.d_1 = d_1
    
    def LHS(self, edge, D_u: float_array, D_v: float_array, k: float) -> complex_array:
        d_1 = self.d_1

        return I_udv(edge, D_u, D_v, k) + d_1/(1j*k)*I_dudv(edge, D_u, D_v, k)

    def RHS(self, edge, D_v: float_array, k: float) -> complex_array:
        raise NotImplementedError("Not implemented yet")



class UltraWeakFlux:
    r"""
    Computes the flux on a inner facet with respect to the degrees
    of freedom from the same cell, that is:
    
    .. math::
        \int_E \left(\left(\varphi_n(\mathbf{x})+\frac{b}{ik}\nabla\varphi_n(\mathbf{x})\cdot\mathbf{n}\right)\overline{\nabla\psi_m(\mathbf{x})\cdot\mathbf{n}}- \left(\vphantom{\frac{1}{2}}aik\varphi_n(\mathbf{x})+\nabla\varphi_n(\mathbf{x})\cdot\mathbf{n}\right)\overline{\psi_m(\mathbf{x})}\right) \,\mathrm{d}S_\mathbf{x}    

    Parameters
    ----------
    edge : Edge or Arc
        Edge parameters.
    D_u : (float, float) array
        Propatagion direction of the trial function.
    D_v : (float, float) array
        Propatagion direction of the test function.
    k : float
        Wavenumber.
    a : float
        Stabilyzing parameter.
    b : float
        Stabilyzing parameter.

    Returns
    -------
    I : complex
        The integral.


    """

    def __init__(self, a: float, b: float):
        self.a = a 
        self.b = b
    
    def LHS(self, edge, D_u: float_array, D_v: float_array, k: float, sign: SIGN) -> complex_array:
        match sign:
            case SIGN.PP:
                a = self.a
                b = self.b
            case SIGN.PM:
                a = self.a
                b = self.b
            case SIGN.MP:
                a = -self.a
                b = -self.b
            case SIGN.MM:
                a = -self.a
                b = -self.b

        k_n = k
        k_m = k

        I = 1/2*I_udv(edge, D_u, D_v, k) + b /(1j*k)*I_dudv(edge, D_u, D_v, k) - a*1j*k*I_uv(edge, D_u, D_v, k) - 1/2*I_duv(edge, D_u, D_v, k)
        match sign:
            case SIGN.PP:
                I = I
            case SIGN.PM:
                I = -I
            case SIGN.MP:
                I = I
            case SIGN.MM:
                I = -I
        return I

# RIGHT NOW RHS only does plane waves as RHS
class DirichletFlux:
    r"""
    Computes the flux on a Neumann boundary condition, that is:

    .. math::
    
        \int_{E}\left(\varphi_n(\mathbf{x})+\frac{d_{1}}{ik}\nabla \varphi_n(\mathbf{x})\cdot\mathbf{n}\right)\overline{\nabla \psi_m(\mathbf{x})\cdot\mathbf{n}}\,\mathrm{d}S_{\mathbf{x}}

    Parameters
    ----------
    edge : Edge or Arc
        Edge parameters.
    D_u : (float, float) array
        Propatagion direction of the trial function.
    D_v : (float, float) array
        Propatagion direction of the test function.
    k : float
        Wavenumber.
    a : float
        Stabilyzing parameter.

    Returns
    -------
    I : complex
        The integral.
    
    """

    def __init__(self, a: float, data = None):
        self.a = a
        self.data = data
    
    def LHS(self, edge, D_u: float_array, D_v: float_array, k: float) -> complex_array:
        a = self.a
        return -I_duv(edge, D_u, D_v, k) + 1j*k*a*I_uv(edge, D_u, D_v, k)
        
    # def RHS(self, edge, D_v: float_array, k: float) -> complex_array:
    #     d_inc = self.data["d_inc"]
    #     a = self.a
    #     return I_uincdv(edge, d_inc, D_v, k) - 1j*a*k*I_uincv(edge, d_inc, D_v, k)
    
    def RHS(self, edge, D_v: float_array, k: float) -> complex_array:
        plane_wave = self.data
        a = self.a
        return I_pw_dv(edge, plane_wave, D_v, k) - 1j*a*k*I_pw_v(edge, plane_wave, D_v, k)