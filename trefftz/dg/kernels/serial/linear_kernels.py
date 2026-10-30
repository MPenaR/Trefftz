from numpy import dot, exp, sinc, pi
from trefftz.numpy_types import float_array
from trefftz.dg.exact import PlaneWave


def I_uv(segment,  k: float, d_u: float_array, d_v: float_array, n_u: float | complex = 1., n_v: float | complex = 1.) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''
    l = segment["l"]
    M = segment["M"]
    T = segment["T"]

    d_uv = n_u*d_u - n_v*d_v

    return l*exp(1j*k*dot(d_uv, M))*sinc(k*l/(2*pi)*dot(d_uv, T))


def I_duv(segment, k: float, d_u: float_array, d_v: float_array, n_u: float | complex = 1., n_v: float | complex = 1.) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''

    N = segment["N"]
    return 1j*k*dot(n_u*d_u, N)*I_uv(segment, k, d_u, d_v, n_u, n_v)


def I_udv(segment, k: float, d_u: float_array, d_v: float_array, n_u: float | complex = 1., n_v: float | complex = 1.) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''

    N = segment["N"]
    return -1j*k*dot(n_v*d_v, N)*I_uv(segment, k, d_u, d_v, n_u, n_v)


def I_dudv(segment, k: float, d_u: float_array, d_v: float_array, n_u: float | complex = 1., n_v: float | complex = 1.) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''
    N = segment["N"]
    return k**2*dot(n_u*d_u, N)*dot(n_v*d_v, N)*I_uv(segment, k, d_u, d_v, n_u, n_v) 

def I_pw_v(segment_v, plane_wave: PlaneWave, k, d_v: float_array, n_v: float | complex = 1.) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{v}\,\mathrm{d}\ell

    where $u_inc$ is a plane wave and $v$ is a plane wave and $E$ is a segment.'''


    return plane_wave.A*I_uv(segment=segment_v, d_u=plane_wave.d, d_v=d_v, k=k, n_v=n_v)


def I_pw_dv(segment_v, plane_wave: PlaneWave, k, d_v: float_array, n_v: float | complex = 1.) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u_inc$ is a plane wave and $v$ is a plane wave and $E$ is a segment.'''

    N = segment_v["N"]

    return -1j*k*dot(N, n_v*d_v)*I_pw_v(segment_v, plane_wave, k, d_v, n_v)


def I_dpw_v(segment_v, plane_wave: PlaneWave, k, d_v: float_array, n_v: float | complex = 1.) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u_inc$ is a plane wave and $v$ is a plane wave and $E$ is a segment.'''

    N = segment_v["N"]

    return 1j*k*dot(N, plane_wave.d)*I_pw_v(segment_v, plane_wave, k, d_v, n_v)

