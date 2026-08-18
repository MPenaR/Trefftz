r"""
Module for the evaluation of fluxes for a single pair of test and trial functions. It was the initial
implementation of the code and now serves as a test for the vectorized implementation. It is tested
against numerical computations of the fluxes.

All around this page, :math:`\Phi_{n}(\mathbf{x})=\exp(ik\mathbf{x}\cdot\mathbf{d}_\varphi)` and :math:`\Psi_{m}(\mathbf{x})=\exp(ik\mathbf{x}\cdot\mathbf{d}_\psi)` 
are the trial and test functions. :math:`l` is the length of an edge, :math:`\mathbf{M}` its midpoint and :math:`\boldsymbol{\tau}` and :math:`\mathbf{n}` 
its tangent and normal unitary vectors. 
"""

from trefftz.numpy_types import float_array
from enum import Enum
from typing import Protocol

from trefftz.dg.kernels.serial_kernels import I_uv, I_duv, I_udv, I_dudv, I_uincdv, I_uincv


JAC_ANGER_MODES = 80

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


class SerialTransmissionKernel(Protocol):
    def LHS(self, edge, d_u: float_array, d_v: float_array, k: float, sign: SIGN) -> complex:
        ...


class SerialLocalKernel(Protocol):
    def LHS(self, edge, d_u: float_array, d_v: float_array, k: float) -> complex:
        ...

    def RHS(self, edge, d_v: float_array, k: float) -> complex:
        ...


class SerialNonLocalKernel(Protocol):
    def LHS(self, edge_u, edge_v, d_u: float_array, d_v: float_array, k: float) -> complex:
        ...


class NeumannFlux:
    r"""
    Computes the flux on a Neumann boundary condition, that is:

    .. math::
    
        \int_{E}\left(\varphi_n(\mathbf{x})+\frac{d_{1}}{ik}\nabla \varphi_n(\mathbf{x})\cdot\mathbf{n}\right)\overline{\nabla \psi_m(\mathbf{x})\cdot\mathbf{n}}\,\mathrm{d}S_{\mathbf{x}}

    This quantity can be exactly evaluated as:

    .. math::
    
        \boxed{-ikl\left(1+d_{1}\mathbf{d}_{n}\cdot\mathbf{n}\right)\mathbf{d}_{m}\cdot\mathbf{n}e^{ik\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\cdot\mathbf{M}}\mathrm{sinc}\left(\frac{kl}{2\pi}\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\cdot\boldsymbol{\tau}\right)}


    Parameters
    ----------
    edge : Edge or Arc
        Edge parameters.
    d_u : (float, float) array
        Propatagion direction of the trial function.
    d_v : (float, float) array
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
    
    def LHS(self, edge, d_u: float_array, d_v: float_array, k: float) -> complex:
        d_1 = self.d_1
        d_m = d_v
        d_n = d_u

        N = edge["N"]

        return I_udv(edge, d_u, d_v, k) + d_1/(1j*k)*I_dudv(edge, d_u, d_v, k)

    def RHS(self, edge, d_v: float_array, k: float) -> complex:
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
    d_u : (float, float) array
        Propatagion direction of the trial function.
    d_v : (float, float) array
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
    
    def LHS(self, edge, d_u: float_array, d_v: float_array, k: float, sign: SIGN) -> complex:
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

        I = 1/2*I_udv(edge, d_u, d_v, k) + b /(1j*k)*I_dudv(edge, d_u, d_v, k) - a*1j*k*I_uv(edge, d_u, d_v, k) - 1/2*I_duv(edge, d_u, d_v, k)
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

    This quantity can be exactly evaluated as:

    .. math::
    
        \boxed{a}


    Parameters
    ----------
    edge : Edge or Arc
        Edge parameters.
    d_u : (float, float) array
        Propatagion direction of the trial function.
    d_v : (float, float) array
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
    
    def LHS(self, edge, d_u: float_array, d_v: float_array, k: float) -> complex:
        a = self.a
        return -I_duv(edge, d_u, d_v, k) - 1j*k*a*I_uv(edge, d_u, d_v, k)
        
    def RHS(self, edge, d_v: float_array, k: float) -> complex:
        d_inc = self.data["d_inc"]
        a = self.a
        return I_uincdv(edge, d_inc, d_v, k) + 1j*a*k*I_uincv(edge, d_inc, d_v, k)