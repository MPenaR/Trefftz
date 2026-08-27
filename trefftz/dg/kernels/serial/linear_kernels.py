from numpy import dot, exp, sinc, pi
from trefftz.numpy_types import float_array
from trefftz.dg.exact import PlaneWave


def I_uv(segment, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''
    l = segment["l"]
    M = segment["M"]
    T = segment["T"]

    return l*exp(1j*k*dot((d_u - d_v), M))*sinc(k*l/(2*pi)*dot(d_u - d_v, T))


def I_duv(segment, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''

    N = segment["N"]
    return 1j*k*dot(d_u, N)*I_uv(segment, d_u, d_v, k)


def I_udv(segment, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''

    N = segment["N"]
    return -1j*k*dot(d_v, N)*I_uv(segment, d_u, d_v, k)


def I_dudv(segment, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''
    N = segment["N"]
    return k**2*dot(d_u, N)*dot(d_v, N)*I_uv(segment, d_u, d_v, k) 

def I_pw_v(segment_v, plane_wave: PlaneWave, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{v}\,\mathrm{d}\ell

    where $u_inc$ is a plane wave and $v$ is a plane wave and $E$ is a segment.'''


    return plane_wave.A*I_uv(segment=segment_v, d_u=plane_wave.d, d_v=d_v, k=k)


def I_pw_dv(segment_v, plane_wave: PlaneWave, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u_inc$ is a plane wave and $v$ is a plane wave and $E$ is a segment.'''

    N = segment_v["N"]

    return -1j*k*dot(N, d_v)*I_pw_v(segment_v, plane_wave, d_v, k)


def I_dpw_v(segment_v, plane_wave: PlaneWave, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u_inc$ is a plane wave and $v$ is a plane wave and $E$ is a segment.'''

    N = segment_v["N"]

    return 1j*k*dot(N, plane_wave.d)*I_pw_v(segment_v, plane_wave, d_v, k)

