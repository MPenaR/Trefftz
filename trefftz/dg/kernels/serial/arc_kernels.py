from numpy import atan2, sin, cos
from numpy.linalg import norm
from trefftz.numpy_types import float_array
from scipy.special import jv
from trefftz.dg.exact import PlaneWave

JAC_ANGER_MODES = 80


def I_uv(arc, d_u: float_array, d_v: float_array, k: float) -> complex:
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


def I_duv(arc, d_u: float_array, d_v: float_array, k: float) -> complex:
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


def I_udv(arc, d_u: float_array, d_v: float_array, k: float) -> complex:
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


def I_dudv(arc, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is an arc of circunference.'''

    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    D_uv = norm(d_u - d_v)
    phi_uv = atan2((d_u - d_v)[1], (d_u - d_v)[0])

    phi_u = atan2(d_u[1], d_u[0])
    phi_v = atan2(d_v[1], d_v[0])

    raise NotImplementedError("need to formally compute this case.")


# def I_uincv(arc, d_inc: float_array, d_v: float_array, k: float) -> complex:
#     r'''Computes the integral:
#     .. math ::
#         \int_E u_{\mathrm{inc}} \overline{v}\,\mathrm{d}\ell

#     where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is an arc of circunference.'''


#     theta_1 = arc["theta_1"]
#     theta_2 = arc["theta_2"]
#     R = arc["R"]

#     D_iv = norm(d_inc - d_v)
#     phi_iv = atan2((d_inc - d_v)[1], (d_inc - d_v)[0])

#     return R*(jv(0, k*R*D_iv)*(theta_2-theta_1) +
#               2*sum(1j**t/t*jv(t, k*R*D_iv)*(sin(t*(theta_2 - phi_iv)) - sin(t*(theta_1 - phi_iv)))
#                     for t in range(1, JAC_ANGER_MODES)))


# def I_uincdv(arc, d_inc: float_array, d_v: float_array, k: float) -> complex:
#     r'''Computes the integral:
#     .. math ::
#         \int_E u_{\mathrm{inc}} \overline{nabla(v)\cdot\mathbf{n}}\,\mathrm{d}\ell

#     where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is an arc of circunference.'''

#     theta_1 = arc["theta_1"]
#     theta_2 = arc["theta_2"]
#     R = arc["R"]

#     D_iv = norm(d_inc - d_v)
#     phi_iv = atan2((d_inc - d_v)[1], (d_inc - d_v)[0])

#     phi_v = atan2(d_v[1], d_v[0])

#     def I(theta: float) -> complex:
#         return -k*R*(-jv(1, k*R*D_iv)*cos(phi_iv - phi_v)*theta + sum(1j**p/p*(jv(p-1, k*R*D_iv)*sin(p*(theta-phi_iv)+(phi_iv-phi_v))-
#                                                                                jv(p+1, k*R*D_iv)*sin(p*(theta-phi_iv)-(phi_iv-phi_v)))
#                                                                                               for p in range(1, JAC_ANGER_MODES)))

#     return I(theta_2) - I(theta_1)


def I_pw_v(arc, plane_wave: PlaneWave, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{v}\,\mathrm{d}\ell
    
    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is an arc of circunference.'''
    
    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    d_iv = plane_wave.d - d_v
    D_iv = norm(d_iv)
    phi_iv = atan2((d_iv)[1], (d_iv)[0])

    I = R*(jv(0, k*R*D_iv)*(theta_2-theta_1) + 2*sum(1j**t/t*jv(t, k*R*D_iv)*(sin(t*(theta_2 - phi_iv)) - sin(t*(theta_1 - phi_iv)))
                                                     for t in range(1, JAC_ANGER_MODES)))
    return plane_wave.A*I

def I_pw_dv(arc, plane_wave: PlaneWave, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{nabla(v)\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is an arc of circunference.'''

    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]
 
    d_iv = plane_wave.d - d_v
    D_iv = norm(d_iv)
    phi_iv = atan2(d_iv[1], d_iv[0])

    phi_v = atan2(d_v[1], d_v[0])

    def I(theta: float) -> complex:
        return -k*R*(-jv(1, k*R*D_iv)*cos(phi_iv - phi_v)*theta + sum(1j**p/p*(jv(p-1, k*R*D_iv)*sin(p*(theta-phi_iv)+(phi_iv-phi_v))-
                                                                               jv(p+1, k*R*D_iv)*sin(p*(theta-phi_iv)-(phi_iv-phi_v)))
                                                                                              for p in range(1, JAC_ANGER_MODES)))

    return plane_wave.A*(I(theta_2) - I(theta_1))
