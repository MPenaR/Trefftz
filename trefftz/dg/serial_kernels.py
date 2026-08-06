from numpy import dot, exp, sinc, pi, atan2, sin, cos
from numpy.linalg import norm
from trefftz.numpy_types import float_array
from typing import Protocol
from scipy.special import jv

JAC_ANGER_MODES = 80

def I_uv(edge, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{v}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is a straight edge.'''
    l = edge["l"]
    M = edge["M"]
    T = edge["T"]

    return l*exp(1j*k*dot((d_u - d_v), M))*sinc(k*l/(2*pi)*dot(d_u - d_v, T))


def I_uv_arc(edge, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{v}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is an arc of circunference.'''

    # to be replaced with edge["theta_2"], edge["theta_1"], and edge["R"] as
    # edge should be a arc_edge

    # P = edge["P"]
    # Q = edge["Q"]

    # theta_1 = atan2(P[1], P[0])
    # theta_2 = atan2(Q[1], Q[0])

    theta_1 = edge["theta_1"]
    theta_2 = edge["theta_2"]
    R = edge["R"]    
    
    D_uv = norm(d_u - d_v)
    phi_uv = atan2((d_u - d_v)[1], (d_u - d_v)[0])

    return R*(jv(0, k*R*D_uv)*(theta_2-theta_1) +
              2*sum(1j**t/t*jv(t, k*R*D_uv)*(sin(t*(theta_2 - phi_uv)) - sin(t*(theta_1 - phi_uv)))
                    for t in range(1, JAC_ANGER_MODES)))

def I_duv_arc(edge, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{v}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is an arc of circunference.'''

    # to be replaced with edge["theta_2"], edge["theta_1"], and edge["R"] as
    # edge should be a arc_edge

    # P = edge["P"]
    # Q = edge["Q"]

    # theta_1 = atan2(P[1], P[0])
    # theta_2 = atan2(Q[1], Q[0])

    theta_1 = edge["theta_1"]
    theta_2 = edge["theta_2"]
    R = edge["R"]

    D_uv = norm(d_u - d_v)
    phi_uv = atan2((d_u - d_v)[1], (d_u - d_v)[0])

    phi_u = atan2(d_u[1], d_u[0])

    primitive = lambda theta: k*R*(-jv(1, k*R*D_uv)*cos(phi_uv - phi_u)* theta + sum(1j**p/p*(jv(p-1,k*R*D_uv)*sin(p*(theta-phi_uv)+(phi_uv-phi_u))-
                                                                                              jv(p+1,k*R*D_uv)*sin(p*(theta-phi_uv)-(phi_uv-phi_u)))
                                                                                              for p in range(1, JAC_ANGER_MODES)))

    return primitive(theta_2) - primitive(theta_1)

# class SerialTransmissionKernel(Protocol):
#     def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float, sign: SIGN) -> complex:
#         ...


# class SerialLocalKernel(Protocol):
#     def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
#         ...

#     def RHS(self, edge, d_psi: float_array, k: float) -> complex:
#         ...


# class SerialNonLocalKernel(Protocol):
#     def LHS(self, edge_u, edge_v, d_phi: float_array, d_psi: float_array, k: float) -> complex:
#         ...

def I_v_arc(edge, d_inc: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{v}\,\mathrm{d}\ell
    
    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is an arc of circunference.'''

    # to be replaced with edge["theta_2"], edge["theta_1"], and edge["R"] as
    # edge should be a arc_edge

    # P = edge["P"]
    # Q = edge["Q"]

    # theta_1 = atan2(P[1], P[0])
    # theta_2 = atan2(Q[1], Q[0])
    
    theta_1 = edge["theta_1"]
    theta_2 = edge["theta_2"]
    R = edge["R"]

    D_iv = norm(d_inc - d_v)
    phi_iv = atan2((d_inc - d_v)[1], (d_inc - d_v)[0])

    return R*(jv(0, k*R*D_iv)*(theta_2-theta_1) +
              2*sum(1j**t/t*jv(t, k*R*D_iv)*(sin(t*(theta_2 - phi_iv)) - sin(t*(theta_1 - phi_iv)))
                    for t in range(1, JAC_ANGER_MODES)))


def I_dv_arc(edge, d_inc: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{nabla(v)\cdot\mathbf{n}}\,\mathrm{d}\ell
    
    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is an arc of circunference.'''

    # to be replaced with edge["theta_2"], edge["theta_1"], and edge["R"] as
    # edge should be a arc_edge

    # P = edge["P"]
    # Q = edge["Q"]

    # theta_1 = atan2(P[1], P[0])
    # theta_2 = atan2(Q[1], Q[0])
    
    theta_1 = edge["theta_1"]
    theta_2 = edge["theta_2"]
    R = edge["R"]

    D_iv = norm(d_inc - d_v)
    phi_iv = atan2((d_inc - d_v)[1], (d_inc - d_v)[0])

    phi_v = atan2(d_v[1], d_v[0])

    primitive = lambda theta: -k*R*(-jv(1, k*R*D_iv)*cos(phi_iv - phi_v)* theta + sum(1j**p/p*(jv(p-1, k*R*D_iv)*sin(p*(theta-phi_iv)+(phi_iv-phi_v))-
                                                                                               jv(p+1, k*R*D_iv)*sin(p*(theta-phi_iv)-(phi_iv-phi_v)))
                                                                                              for p in range(1, JAC_ANGER_MODES)))

    return primitive(theta_2) - primitive(theta_1)
