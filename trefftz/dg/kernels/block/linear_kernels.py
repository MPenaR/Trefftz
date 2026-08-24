from numpy import dot, exp, sinc, pi
from trefftz.numpy_types import float_array, complex_array
from trefftz.dg.exact import PlaneWave

def I_uv(segment, D_u: float_array, D_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''
    l = segment["l"]
    M = segment["M"]
    T = segment["T"]

    d_uv = D_u[None, : , : ] - D_v[ :, None, : ]

    return l*exp(1j*k*dot(d_uv, M))*sinc(k*l/(2*pi)*dot(d_uv, T))


def I_duv(segment, D_u: float_array, D_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''

    N = segment["N"]

    return 1j*k*dot(D_u, N)[None, :]*I_uv(segment, D_u, D_v, k)


def I_udv(segment, D_u: float_array, D_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''

    N = segment["N"]
    return -1j*k*dot(D_v, N)[:, None]*I_uv(segment, D_u, D_v, k)


def I_dudv(segment, D_u: float_array, D_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''
    N = segment["N"]
    return k**2*dot(D_u, N)[None, :]*dot(D_v, N)[:, None]*I_uv(segment, D_u, D_v, k)


# def I_uincv(segment, d_inc: float_array, D_v: float_array, k: float) -> complex_array:
#     r'''Computes the integral:
#     .. math ::
#         \int_E u_{\mathrm{inc}} \overline{v}\,\mathrm{d}\ell

#     where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is a segment.'''

#     l = segment["l"]
#     M = segment["M"]
#     T = segment["T"]

#     d_incv = d_inc[None, : ] - D_v[ :, : ]

#     return l*exp(1j*k*dot(d_incv, M))*sinc(k*l/(2*pi)*dot(d_incv, T))


# def I_uincdv(segment, d_inc: float_array, D_v: float_array, k: float) -> complex_array:
#     r'''Computes the integral:
#     .. math ::
#         \int_E u_{\mathrm{inc}} \overline{\nabla v \cdot \mathbf{n}}\,\mathrm{d}\ell

#     where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is a segment.'''

#     N = segment["N"]
#     return -1j*k*dot(D_v, N)*I_uincv(segment, d_inc, D_v, k)


def I_pw_v(segment_v, plane_wave: PlaneWave, D_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{v}\,\mathrm{d}\ell

    where $u_inc$ is a plane wave and $v$ is a plane wave and $E$ is a segment.'''
    l = segment_v["l"]
    M = segment_v["M"]
    T = segment_v["T"]

    d_iv = plane_wave.d - D_v

    return plane_wave.A*l*exp(1j*k*dot(d_iv, M))*sinc(k*l/(2*pi)*dot(d_iv, T))


def I_pw_dv(segment_v, plane_wave: PlaneWave, D_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u_inc$ is a plane wave and $v$ is a plane wave and $E$ is a segment.'''

    N = segment_v["N"]

    return -1j*k*dot(D_v, N)*I_pw_v(segment_v, plane_wave, D_v, k)
