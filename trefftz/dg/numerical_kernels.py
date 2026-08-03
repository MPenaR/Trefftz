from numpy import linspace, outer, sin, cos, pi, exp, dot, conj, isclose, array
from numpy.lib.scimath import sqrt
from numpy.linalg import norm
from numpy import trapezoid as Int
import numpy as np
from scipy.special import hankel1, h1vp
from trefftz.numpy_types import float_array, complex_array

N_POINTS = int(1E5)
NtD_MODES = 20


class numNeumannKernel:

    r"""Implements the integral: 
    .. math ::

        \int_{E}\left(u+\frac{\mathfrak{d}_{1}}{ik}\nabla u\cdot\mathbf{n}\right)\overline{\nabla v\cdot n}d\ell
    
    """
    def __init__(self, d_1: float):
        self.d_1 = d_1

    
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        d1 = self.d_1
        d_n = d_phi
        d_m = d_psi

        t = linspace(0, 1, N_POINTS)
        P = edge["P"]
        T = edge["T"]
        l = edge["l"]
        N = edge["N"]
        x = P + outer(t, T)*l
        u = exp(1j*k*dot(x, d_n))
        v = exp(1j*k*dot(x, d_m))
        du_dn = 1j*k*dot(N, d_n)*u
        dv_dn = 1j*k*dot(N, d_m)*v

        I = l*Int((u + d1/(1j*k)*du_dn)*conj(dv_dn), t)
        return I

    def RHS(self) -> complex:
        raise NotImplementedError
        
class numUltraWeakKernel:

    r"""Implements the integral: 
    .. math ::

        \int_{E}\left(\left(u +\frac{\mathfrak{b}}{ik}\nabla u\cdot\mathbf{n}\right)\overline{\nabla v\cdot\mathbf{n}}-\left(\mathfrak{a}iku +\nabla u \cdot \mathbf{n}\right)\overline{v}\right)d\ell
    
    """
    def __init__(self, a: float, b: float):
        self.a = a
        self.b = b

    
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        a = self.a
        b = self.b

        d_n = d_phi
        d_m = d_psi

        t = linspace(0, 1, N_POINTS)
        P = edge["P"]
        T = edge["T"]
        l = edge["l"]
        N = edge["N"]
        x = P + outer(t, T)*l
        u = exp(1j*k*dot(x, d_n))
        v = exp(1j*k*dot(x, d_m))
        du_dn = 1j*k*dot(N, d_n)*u
        dv_dn = 1j*k*dot(N, d_m)*v

        I = l*Int((u/2 + b/(1j*k)*du_dn)*conj(dv_dn) - (a*1j*k*u + du_dn/2)*conj(v), t)
        return I


def numI_uv_arc(edge, d_u: float_array, d_v: float_array, k: float, R: float) -> complex:
    P = edge["P"]
    Q = edge["Q"]
    theta_1 = np.atan2(P[1], P[0])
    theta_2 = np.atan2(Q[1], Q[0])

    theta = np.linspace(theta_1, theta_2, N_POINTS)
    u_r = np.column_stack([np.cos(theta), np.sin(theta)])
    x = R*u_r
    u = exp(1j*k*dot(x, d_u))
    v = exp(1j*k*dot(x, d_v))
    I = Int(u*conj(v), theta)*R  # int u*conj(v) dl 
    return I


def numI_uv(edge, d_u: float_array, d_v: float_array, k: float) -> complex:
    P = edge["P"]
    T = edge["T"]
    l = edge["l"]
    t = np.linspace(0, 1, N_POINTS)

    x = P + l*outer(t,T)
    u = exp(1j*k*dot(x, d_u))
    v = exp(1j*k*dot(x, d_v))
    I = Int(u*conj(v), t)*l  # int u*conj(v) dl 
    return I


def numI_duv_arc(edge, d_u: float_array, d_v: float_array, k: float, R: float) -> complex:
    P = edge["P"]
    Q = edge["Q"]
    theta_1 = np.atan2(P[1], P[0])
    theta_2 = np.atan2(Q[1], Q[0])

    theta = np.linspace(theta_1, theta_2, N_POINTS)
    u_r = np.column_stack([np.cos(theta), np.sin(theta)])
    N = u_r
    x = R*u_r

    u = exp(1j*k*dot(x, d_u))
    v = exp(1j*k*dot(x, d_v))

    du_dn = 1j*k*dot(N, d_u)*u

    I = Int(du_dn*conj(v), theta)*R  # int du/dn*conj(v) dl 
    return I
