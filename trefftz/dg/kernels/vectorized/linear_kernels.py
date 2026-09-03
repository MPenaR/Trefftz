from numpy import dot, exp, sinc, pi, einsum
from trefftz.numpy_types import float_array, complex_array
from trefftz.dg.exact import PlaneWave

def I_uv(segments, D_u: float_array, n_u: float | complex,  D_v: float_array, n_v: float | complex,  k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ are several edges at once.'''
    l = segments["l"] # (N_E)
    M = segments["M"] # (N_E, 2)
    T = segments["T"] # (N_E, 2)

    # D_u: (N_E, N_theta, 2)
    # D_v: (N_E, N_theta, 2)
    # n_u: (N_E)
    # n_v: (N_E)

    d_uv = n_u[:, None, None , None ]*D_u[None, None, : , : ] - n_v[:, None, None , None ]*D_v[ None, :, None, : ] # (N_E, N_theta, N_theta, 2)

    # return l*exp(1j*k*dot(d_uv, M))*sinc(k*l/(2*pi)*dot(d_uv, T))   # (N_E, N_theta, N_theta)
    return l[:, None, None]*exp(1j*k*einsum("eijt,et->eij", d_uv, M))*sinc(k*l[:, None, None]/(2*pi)*einsum("eijt,et->eij", d_uv, T))   # (N_E, N_theta, N_theta)


def I_duv(segment, D_u: float_array, n_u: float | complex,  D_v: float_array, n_v: float | complex, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''

    N = segment["N"]

    return 1j*k*dot(n_u*D_u, N)[None, :]*I_uv(segment, D_u, n_u, D_v, n_v, k)


def I_udv(segment, D_u: float_array, n_u: float | complex,  D_v: float_array, n_v: float | complex, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''

    N = segment["N"]
    return -1j*k*dot(n_v*D_v, N)[:, None]*I_uv(segment, D_u, n_u, D_v, n_v, k)


def I_dudv(segment, D_u: float_array, n_u: float | complex,  D_v: float_array, n_v: float | complex, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E nabla(u)\cdot\mathbf{n} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is a segment.'''
    N = segment["N"]
    return k**2*dot(n_u*D_u, N)[None, :]*dot(n_v*D_v, N)[:, None]*I_uv(segment, D_u, n_u, D_v, n_v, k)


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
