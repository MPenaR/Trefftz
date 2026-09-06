from numpy import atan2, sin, cos
from numpy.linalg import norm
from trefftz.numpy_types import float_array, complex_array
from scipy.special import jv
from trefftz.dg.exact import PlaneWave

JAC_ANGER_MODES = 80


def I_uv(arcs, D_u: float_array, n_u: float_array, D_v: float_array, n_v: float_array,  k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ are a set of arcs of circunference.'''

    theta_1 = arcs["theta_1"]
    theta_2 = arcs["theta_2"]
    R = arcs["R"]


    # n_u: (N_E)
    # n_v: (N_E)

    d_uv = n_u[:, None, None, None]*D_u[None, None, :, :] - n_v[:, None, None, None]*D_v[None, :, None, :]
    D_uv = norm(d_uv, axis=-1)
    phi_uv = atan2(d_uv[:,  :, :, 1], d_uv[:, :, :, 0])

    return R[:, None, None]*(jv(0, k*R[:, None, None]*D_uv)*(theta_2-theta_1)[:, None, None] +
              2*sum(1j**t/t*jv(t, k*R[:, None, None]*D_uv)*(sin(t*(theta_2[:, None, None] - phi_uv)) - sin(t*(theta_1[:, None, None] - phi_uv)))
                    for t in range(1, JAC_ANGER_MODES)))


def I_duv(arcs, D_u: float_array, n_u: float_array, D_v: float_array, n_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is an arc of circunference.'''

    theta_1 = arcs["theta_1"]
    theta_2 = arcs["theta_2"]
    R = arcs["R"]

    d_uv = n_u[:, None, None, None]*D_u[None, None, :, :] - n_v[:, None, None, None]*D_v[None, :, None, :]
    D_uv = norm(d_uv, axis=-1)
    phi_uv = atan2(d_uv[:,  :, :, 1], d_uv[:, :, :, 0])

    phi_u = atan2(D_u[ :, 1], D_u[ :, 0])

    def I(theta: float_array) -> complex_array:
        return k*n_u[:, None, None]*R[:, None, None]*(-jv(1, k*R[:, None, None]*D_uv)*cos(phi_uv - phi_u[None, None, :])*theta[:, None, None] + sum(1j**p/p*(jv(p-1, k*R[:, None, None]*D_uv)*sin(p*(theta[:, None, None]-phi_uv)+(phi_uv-phi_u[None, None, :]))-
                                                                                                                                                             jv(p+1, k*R[:, None, None]*D_uv)*sin(p*(theta[:, None, None]-phi_uv)-(phi_uv-phi_u[None, None, :])))
                                                                                              for p in range(1, JAC_ANGER_MODES))) 
    return I(theta_2) - I(theta_1)


def I_udv(arcs, D_u: float_array, n_u: float_array, D_v: float_array, n_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is an arc of circunference.'''

    theta_1 = arcs["theta_1"]
    theta_2 = arcs["theta_2"]
    R = arcs["R"]

    d_uv = n_u[:, None, None, None]*D_u[None, None, :, :] - n_v[:, None, None, None]*D_v[None, :, None, :]
    D_uv = norm(d_uv, axis=-1)
    phi_uv = atan2(d_uv[:,  :, :, 1], d_uv[:, :, :, 0])

    phi_v = atan2(D_v[:, 1], D_v[:, 0])

    def I(theta: float) -> complex_array:
        return -k*n_v[:, None, None]*R[:, None, None]*(-jv(1, k*R[:, None, None]*D_uv)*cos(phi_uv - phi_v[None, :, None])* theta[:, None, None] + sum(1j**p/p*(jv(p-1, k*R[:, None, None]*D_uv)*sin(p*(theta[:, None, None]-phi_uv)+(phi_uv-phi_v[None, :, None]))-
                                                                                                                                                jv(p+1, k*R[:, None, None]*D_uv)*sin(p*(theta[:, None, None]-phi_uv)-(phi_uv-phi_v[None, :, None])))
                                                                                              for p in range(1, JAC_ANGER_MODES)))

    return I(theta_2) - I(theta_1)


def I_dudv(arcs, D_u: float_array, n_u: float | complex, D_v: float_array, n_v: float | complex, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is an arc of circunference.'''

    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    d_uv = n_u*D_u[None, :, :] - n_v*D_v[:, None, :]
    D_uv = norm(d_uv, axis=-1)
    phi_uv = atan2(d_uv[ :, :, 1], d_uv[ :, :, 0])

    phi_u = atan2(D_v[:, 1], D_u[:, 0])
    phi_v = atan2(D_v[:, 1], D_v[:, 0])

    raise NotImplementedError("waiting for the serial version")


def I_pw_v(arc, plane_wave: PlaneWave, D_v: float_array, n_v: float | complex, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{v}\,\mathrm{d}\ell
    
    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is an arc of circunference.'''
    
    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]
     
    d_iv = plane_wave.d - n_v*D_v 
    D_iv = norm(d_iv, axis=-1)
    phi_iv = atan2(d_iv[:, 1], d_iv[:, 0])

    I = R*(jv(0, k*R*D_iv)*(theta_2-theta_1) +
              2*sum(1j**t/t*jv(t, k*R*D_iv)*(sin(t*(theta_2 - phi_iv)) - sin(t*(theta_1 - phi_iv)))
                    for t in range(1, JAC_ANGER_MODES)))
    return plane_wave.A*I

def I_pw_dv(arc, plane_wave: PlaneWave, D_v: float_array, n_v: float | complex, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{nabla(v)\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is an arc of circunference.'''

    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    d_iv = plane_wave.d - n_v*D_v 
    D_iv = norm(d_iv, axis=-1)
    phi_iv = atan2(d_iv[:, 1], d_iv[:, 0])

    phi_v = atan2(D_v[:, 1], D_v[:, 0])

    def I(theta: float) -> complex_array:
        return -k*n_v*R*(-jv(1, k*R*D_iv)*cos(phi_iv - phi_v)*theta + sum(1j**p/p*(jv(p-1, k*R*D_iv)*sin(p*(theta-phi_iv)+(phi_iv-phi_v))-
                                                                               jv(p+1, k*R*D_iv)*sin(p*(theta-phi_iv)-(phi_iv-phi_v)))
                                                                                              for p in range(1, JAC_ANGER_MODES)))

    return plane_wave.A*(I(theta_2) - I(theta_1))


def I_dpw_v(arc, plane_wave: PlaneWave, D_v: float_array, n_v: float | complex, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is an arc of circunference.'''

    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    d_iv = plane_wave.d[None, :] - n_v*D_v
    D_iv = norm(d_iv, axis=-1)
    phi_iv = atan2(d_iv[ :, 1], d_iv[ :, 0])

    phi_i = atan2(plane_wave.d[1], plane_wave.d[0])

    def I(theta: float) -> complex_array:
        return k*R*(-jv(1, k*R*D_iv)*cos(phi_iv - phi_i)*theta + sum(1j**p/p*(jv(p-1, k*R*D_iv)*sin(p*(theta-phi_iv)+(phi_iv-phi_i))-
                                                                                       jv(p+1, k*R*D_iv)*sin(p*(theta-phi_iv)-(phi_iv-phi_i)))
                                                                                              for p in range(1, JAC_ANGER_MODES))) 
    return I(theta_2) - I(theta_1)
