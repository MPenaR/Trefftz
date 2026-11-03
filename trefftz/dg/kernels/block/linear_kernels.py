from numpy import dot, exp, sinc, pi
from trefftz.numpy_types import float_array, complex_array
from trefftz.dg.exact import PlaneWave

def I_uv(segment,  k: float, D_u: float_array,  D_v: float_array, n_u: float | complex = 1., n_v: float | complex = 1.) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''
    l = segment["l"]
    M = segment["M"]
    T = segment["T"]

    d_uv = n_u*D_u[None, : , : ] - n_v*D_v[ :, None, : ]

    return l*exp(1j*k*dot(d_uv, M))*sinc(k*l/(2*pi)*dot(d_uv, T))


def I_duv(segment,  k: float, D_u: float_array,  D_v: float_array, n_u: float | complex = 1., n_v: float | complex = 1.) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''

    N = segment["N"]

    return 1j*k*dot(n_u*D_u, N)[None, :]*I_uv(segment, k, D_u, D_v, n_u, n_v)


def I_udv(segment,  k: float, D_u: float_array,  D_v: float_array, n_u: float | complex = 1., n_v: float | complex = 1.) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''

    N = segment["N"]
    return -1j*k*dot(n_v*D_v, N)[:, None]*I_uv(segment, k, D_u, D_v, n_u, n_v)


def I_dudv(segment,  k: float, D_u: float_array,  D_v: float_array, n_u: float | complex = 1., n_v: float | complex = 1.) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''
    N = segment["N"]
    return k**2*dot(n_u*D_u, N)[None, :]*dot(n_v*D_v, N)[:, None]*I_uv(segment, k, D_u, D_v, n_u, n_v)


def I_pw_v(segment_v, plane_wave: PlaneWave, D_v: float_array, n_v: float | complex,  k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{v}\,\mathrm{d}\ell

    where $u_inc$ is a plane wave and $v$ is a plane wave and $E$ is a segment.'''
    l = segment_v["l"]
    M = segment_v["M"]
    T = segment_v["T"]

    d_iv = plane_wave.d - n_v*D_v

    return plane_wave.A*l*exp(1j*k*dot(d_iv, M))*sinc(k*l/(2*pi)*dot(d_iv, T))


def I_pw_dv(segment_v, plane_wave: PlaneWave, D_v: float_array, n_v: float | complex, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u_inc$ is a plane wave and $v$ is a plane wave and $E$ is a segment.'''

    N = segment_v["N"]

    return -1j*k*dot(n_v*D_v, N)*I_pw_v(segment_v, plane_wave, D_v, n_v, k)


def I_dpw_v(segment_v, plane_wave: PlaneWave, D_v: float_array, n_v: float | complex, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u_inc$ is a plane wave and $v$ is a plane wave and $E$ is a segment.'''

    N = segment_v["N"]

    return 1j*k*dot(plane_wave.d, N)*I_pw_v(segment_v, plane_wave, D_v, n_v, k)
