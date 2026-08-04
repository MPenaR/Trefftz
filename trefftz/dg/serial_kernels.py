from numpy import dot, exp, sinc, pi, atan2, sin, cos
from numpy.linalg import norm
from trefftz.numpy_types import float_array
from enum import Enum
from typing import Protocol
from scipy.special import jv

JAC_ANGER_MODES = 30

def I_uv(edge, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{v}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is a straight edge.'''
    l = edge["l"]
    M = edge["M"]
    T = edge["T"]

    return l*exp(1j*k*dot((d_u - d_v), M))*sinc(k*l/(2*pi)*dot(d_u - d_v, T))


def I_uv_arc(edge, d_u: float_array, d_v: float_array, k: float, R: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{v}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is an arc of circunference.'''

    # to be replaced with edge["theta_2"], edge["theta_1"], and edge["R"] as
    # edge should be a arc_edge

    P = edge["P"]
    Q = edge["Q"]

    theta_1 = atan2(P[1], P[0])
    theta_2 = atan2(Q[1], Q[0])
    
    D_uv = norm(d_u - d_v)
    phi_uv = atan2((d_u - d_v)[1], (d_u - d_v)[0])

    return R*(jv(0, k*R*D_uv)*(theta_2-theta_1) +
              2*sum(1j**t/t*jv(t, k*R*D_uv)*(sin(t*(theta_2 - phi_uv)) - sin(t*(theta_1 - phi_uv)))
                    for t in range(1, JAC_ANGER_MODES)))

def I_duv_arc(edge, d_u: float_array, d_v: float_array, k: float, R: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{v}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is an arc of circunference.'''

    # to be replaced with edge["theta_2"], edge["theta_1"], and edge["R"] as
    # edge should be a arc_edge

    P = edge["P"]
    Q = edge["Q"]

    theta_1 = atan2(P[1], P[0])
    theta_2 = atan2(Q[1], Q[0])
    
    D_uv = norm(d_u - d_v)
    phi_uv = atan2((d_u - d_v)[1], (d_u - d_v)[0])

    phi_u = atan2(d_u[1], d_u[0])

    primitive = lambda theta: k*R*(-jv(1, k*R*D_uv)*cos(phi_uv - phi_u)* theta + sum(1j**p/p*(jv(p-1,k*R*D_uv)*sin(p*(theta-phi_uv)+(phi_uv-phi_u))-
                                                                                              jv(p+1,k*R*D_uv)*sin(p*(theta-phi_uv)-(phi_uv-phi_u)))
                                                                                              for p in range(1, JAC_ANGER_MODES)))

    return primitive(theta_2) - primitive(theta_1)

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
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float, sign: SIGN) -> complex:
        ...


class SerialLocalKernel(Protocol):
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        ...

    def RHS(self, edge, d_psi: float_array, k: float) -> complex:
        ...


class SerialNonLocalKernel(Protocol):
    def LHS(self, edge_u, edge_v, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        ...


class NeumannKernel:
    '''Serial Neumann kernel'''
    def __init__(self, d_1: float):
        self.d_1 = d_1
    
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        d_1 = self.d_1
        d_m = d_psi
        d_n = d_phi

        N = edge["N"]

        return -1j*k*(1 + d_1 * dot(d_n, N))*dot(d_m, N)*I_uv(edge, d_phi, d_psi, k)


    def RHS(self, edge, d_psi: float_array, k: float) -> complex:
        raise NotImplementedError("Not implemented yet")


class DirichletKernel:
    '''Serial Dirichlet kernel'''
    def __init__(self, a: float, data = None):
        self.a = a
        self.data = data
    
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        a = self.a
        d_n = d_phi
        N = edge["N"]

        return -1j*k*(dot(d_n, N) + a)*I_uv(edge, d_phi, d_psi, k)
        
    def RHS(self, edge, d_psi: float_array, k: float) -> complex:
        raise NotImplementedError("Not implemented yet")


class CircularDirichletKernel:
    '''Serial Dirichlet kernel'''
    def __init__(self, a: float, R: float, data = None):
        self.a = a
        self.R = R
        self.data = data
    
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        a = self.a
        R = self.R

        return -I_duv_arc(edge, d_phi, d_psi, k, R) - 1j*k*a*I_uv_arc(edge, d_phi, d_psi, k, R)
        
    def RHS(self, edge, d_psi: float_array, k: float) -> complex:
        d_inc = self.data["d_inc"]
        a = self.a
        R = self.R
        return I_dv_arc(edge, d_inc, d_psi, k, R) + 1j*a*k*I_v_arc(edge, d_inc, d_psi, k, R)





class UltraWeakKernel:
    '''Transmission kernel for the UWVF'''
    def __init__(self, a: float, b: float):
        self.a = a 
        self.b = b
    
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float, sign: SIGN) -> complex:
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

        d_m = d_psi
        d_n = d_phi

        k_n = k
        k_m = k

        M = edge["M"]
        l = edge["l"]
        N = edge["N"]
        T = edge["T"]

        # I = -1j*k*l*((1/2+b*dot(d_n, N))*dot(d_m, N) + 1/2*dot(d_n, N) + a)*exp(1j*k*dot(d_n-d_m, M))*sinc(k*l/(2*pi)*dot(d_n-d_m, T))
        I = -1j*k*((1/2+b*dot(d_n, N))*dot(d_m, N) + 1/2*dot(d_n, N) + a)*I_uv(edge, d_n, d_m, k)

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
    

def I_v_arc(edge, d_inc: float_array, d_v: float_array, k: float, R: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{v}\,\mathrm{d}\ell
    
    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is an arc of circunference.'''

    # to be replaced with edge["theta_2"], edge["theta_1"], and edge["R"] as
    # edge should be a arc_edge

    P = edge["P"]
    Q = edge["Q"]

    theta_1 = atan2(P[1], P[0])
    theta_2 = atan2(Q[1], Q[0])
    
    D_iv = norm(d_inc - d_v)
    phi_iv = atan2((d_inc - d_v)[1], (d_inc - d_v)[0])

    return R*(jv(0, k*R*D_iv)*(theta_2-theta_1) +
              2*sum(1j**t/t*jv(t, k*R*D_iv)*(sin(t*(theta_2 - phi_iv)) - sin(t*(theta_1 - phi_iv)))
                    for t in range(1, JAC_ANGER_MODES)))

def I_dv_arc(edge, d_inc: float_array, d_v: float_array, k: float, R: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{nabla(v)\cdot\mathbf{n}}\,\mathrm{d}\ell
    
    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is an arc of circunference.'''

    # to be replaced with edge["theta_2"], edge["theta_1"], and edge["R"] as
    # edge should be a arc_edge

    P = edge["P"]
    Q = edge["Q"]

    theta_1 = atan2(P[1], P[0])
    theta_2 = atan2(Q[1], Q[0])
    
    D_iv = norm(d_inc - d_v)
    phi_iv = atan2((d_inc - d_v)[1], (d_inc - d_v)[0])

    phi_v = atan2(d_v[1], d_v[0])

    primitive = lambda theta: -k*R*(-jv(1, k*R*D_iv)*cos(phi_iv - phi_v)* theta + sum(1j**p/p*(jv(p-1, k*R*D_iv)*sin(p*(theta-phi_iv)+(phi_iv-phi_v))-
                                                                                               jv(p+1, k*R*D_iv)*sin(p*(theta-phi_iv)-(phi_iv-phi_v)))
                                                                                              for p in range(1, JAC_ANGER_MODES)))

    return primitive(theta_2) - primitive(theta_1)
