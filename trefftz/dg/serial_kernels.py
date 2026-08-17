from numpy import dot, exp, sinc, pi, atan2, sin, cos
from numpy.linalg import norm
from trefftz.numpy_types import float_array
from trefftz.mesh.core2 import edge_dtype, arc_dtype
from scipy.special import jv

JAC_ANGER_MODES = 80


def I_uv(edge, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{v}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is either an arc of circunference or a segment.'''
    if edge.dtype == edge_dtype:
        return I_uv_segment(edge, d_u, d_v, k)
    elif edge.dtype == arc_dtype:
        return I_uv_arc(edge, d_u, d_v, k)


def I_duv(edge, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E \nabla u\cdot\mathbf{n} \overline{v}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is either an arc of circunference or a segment.'''
    if edge.dtype == edge_dtype:
        return I_duv_segment(edge, d_u, d_v, k)
    elif edge.dtype == arc_dtype:
        return I_duv_arc(edge, d_u, d_v, k)


def I_udv(edge, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{\nabla v \cdot \mathbf{n}}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is either an arc of circunference or a segment.'''
    if edge.dtype == edge_dtype:
        return I_udv_segment(edge, d_u, d_v, k)
    elif edge.dtype == arc_dtype:
        return I_udv_arc(edge, d_u, d_v, k)


def I_dudv(edge, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E \nabla u \cdot\mathbf{n} \overline{\nabla v \cdot \mathbf{n}}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is either an arc of circunference or a segment.'''
    if edge.dtype == edge_dtype:
        return I_dudv_segment(edge, d_u, d_v, k)
    elif edge.dtype == arc_dtype:
        return I_dudv_arc(edge, d_u, d_v, k)



def I_uincv(edge, d_inc: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
.. math ::
        \int_E u_{\mathrm{inc}} \overline{v}\,\mathrm{d}\ell
    
    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is either an arc of circunference or a segment.'''
    if edge.dtype == edge_dtype:
        return I_uincv_segment(edge, d_inc, d_v, k)
    elif edge.dtype == arc_dtype:
        return I_uincv_arc(edge, d_inc, d_v, k)


def I_uincdv(edge, d_inc: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell
    
    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is either an arc of circunference or a segment.'''
    if edge.dtype == edge_dtype:
        return I_uincdv_segment(edge, d_inc, d_v, k)
    elif edge.dtype == arc_dtype:
        return I_uincdv_arc(edge, d_inc, d_v, k)



def I_uv_segment(segment, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{v}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is a segment.'''
    l = segment["l"]
    M = segment["M"]
    T = segment["T"]

    return l*exp(1j*k*dot((d_u - d_v), M))*sinc(k*l/(2*pi)*dot(d_u - d_v, T))


def I_uv_arc(arc, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{v}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is an arc of circunference.'''


    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]    
    
    D_uv = norm(d_u - d_v)
    phi_uv = atan2((d_u - d_v)[1], (d_u - d_v)[0])

    return R*(jv(0, k*R*D_uv)*(theta_2-theta_1) +
              2*sum(1j**t/t*jv(t, k*R*D_uv)*(sin(t*(theta_2 - phi_uv)) - sin(t*(theta_1 - phi_uv)))
                    for t in range(1, JAC_ANGER_MODES)))


def I_duv_segment(segment, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{v}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is a segment.'''
    N = segment["N"]
    return 1j*k*dot(d_u, N)*I_uv_segment(segment, d_u, d_v, k)


def I_duv_arc(arc, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{v}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is an arc of circunference.'''

    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    D_uv = norm(d_u - d_v)
    phi_uv = atan2((d_u - d_v)[1], (d_u - d_v)[0])

    phi_u = atan2(d_u[1], d_u[0])

    def I(theta: float) -> complex:
        return k*R*(-jv(1, k*R*D_uv)*cos(phi_uv - phi_u)* theta + sum(1j**p/p*(jv(p-1, k*R*D_uv)*sin(p*(theta-phi_uv)+(phi_uv-phi_u))-
                                                                               jv(p+1, k*R*D_uv)*sin(p*(theta-phi_uv)-(phi_uv-phi_u)))
                                                                                              for p in range(1, JAC_ANGER_MODES))) 
    return I(theta_2) - I(theta_1)


def I_udv_segment(segment, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is a segment.'''
    N = segment["N"]
    return -1j*k*dot(d_v, N)*I_uv_segment(segment, d_u, d_v, k)


def I_udv_arc(arc, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is an arc of circunference.'''

    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    D_uv = norm(d_u - d_v)
    phi_uv = atan2((d_u - d_v)[1], (d_u - d_v)[0])

    phi_v = atan2(d_v[1], d_v[0])

    def I(theta: float) -> complex:
        return -k*R*(-jv(1, k*R*D_uv)*cos(phi_uv - phi_v)* theta + sum(1j**p/p*(jv(p-1, k*R*D_uv)*sin(p*(theta-phi_uv)+(phi_uv-phi_v))-
                                                                                jv(p+1, k*R*D_uv)*sin(p*(theta-phi_uv)-(phi_uv-phi_v)))
                                                                                              for p in range(1, JAC_ANGER_MODES)))

    return I(theta_2) - I(theta_1)

def I_dudv_segment(segment, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is a segment.'''
    N = segment["N"]
    return k**2*dot(d_u, N)*dot(d_v, N)*I_uv_segment(segment, d_u, d_v, k)



def I_dudv_arc(arc, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell
    
    where $u$ and $v$ are plane waves and $E$ is an arc of circunference.'''

    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    D_uv = norm(d_u - d_v)
    phi_uv = atan2((d_u - d_v)[1], (d_u - d_v)[0])

    phi_u = atan2(d_u[1], d_u[0])
    phi_v = atan2(d_v[1], d_v[0])

    raise NotImplementedError("not implemented yet, needs int cos cos exp")
    # def I(theta: float) -> complex: 
    #     return ...
    
    # return I(theta_2) - I(theta_1)


def I_uincv_segment(segment, d_inc: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{v}\,\mathrm{d}\ell
    
    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is a segment.'''

    l = segment["l"]
    M = segment["M"]
    T = segment["T"]

    return l*exp(1j*k*dot((d_inc - d_v), M))*sinc(k*l/(2*pi)*dot(d_inc - d_v, T))

def I_uincv_arc(arc, d_inc: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{v}\,\mathrm{d}\ell
    
    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is an arc of circunference.'''

    
    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    D_iv = norm(d_inc - d_v)
    phi_iv = atan2((d_inc - d_v)[1], (d_inc - d_v)[0])

    return R*(jv(0, k*R*D_iv)*(theta_2-theta_1) +
              2*sum(1j**t/t*jv(t, k*R*D_iv)*(sin(t*(theta_2 - phi_iv)) - sin(t*(theta_1 - phi_iv)))
                    for t in range(1, JAC_ANGER_MODES)))


def I_uincdv_segment(segment, d_inc: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{\nabla v \cdot \mathbf{n}}\,\mathrm{d}\ell
    
    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is a segment.'''

    l = segment["l"]
    M = segment["M"]
    T = segment["T"]
    N = segment["N"]
    return -1j*k*dot(d_v, N)*I_uv_segment(segment, d_inc, d_v, k)


def I_uincdv_arc(arc, d_inc: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{nabla(v)\cdot\mathbf{n}}\,\mathrm{d}\ell
    
    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is an arc of circunference.'''

    
    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    D_iv = norm(d_inc - d_v)
    phi_iv = atan2((d_inc - d_v)[1], (d_inc - d_v)[0])

    phi_v = atan2(d_v[1], d_v[0])

    primitive = lambda theta: -k*R*(-jv(1, k*R*D_iv)*cos(phi_iv - phi_v)* theta + sum(1j**p/p*(jv(p-1, k*R*D_iv)*sin(p*(theta-phi_iv)+(phi_iv-phi_v))-
                                                                                               jv(p+1, k*R*D_iv)*sin(p*(theta-phi_iv)-(phi_iv-phi_v)))
                                                                                              for p in range(1, JAC_ANGER_MODES)))

    return primitive(theta_2) - primitive(theta_1)
