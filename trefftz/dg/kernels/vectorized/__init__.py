from trefftz.numpy_types import float_array, complex_array
from trefftz.mesh.core import edge_dtype, arc_dtype
from trefftz.dg.kernels.vectorized import arc_kernels, linear_kernels
from trefftz.dg.exact import PlaneWave
'''module for vectorized evaluation of kernels over arbitrary sets of edges'''


_KERNELS = {
    edge_dtype: linear_kernels,
    arc_dtype: arc_kernels
}


def I_uv(edges, D_u: float_array, n_u: float_array, D_v: float_array, n_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is either an arc of circunference or a segment.'''

    return _KERNELS[edges.dtype].I_uv(edges, D_u, n_u, D_v, n_v, k)



def I_duv(edges, D_u: float_array, n_u: float_array, D_v: float_array, n_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E \nabla u\cdot\mathbf{n} \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is either an arc of circunference or a segment.'''

    return _KERNELS[edges.dtype].I_duv(edges, D_u, n_u, D_v, n_v, k)


def I_udv(edges, D_u: float_array, n_u: float_array, D_v: float_array, n_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{\nabla v \cdot \mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is either an arc of circunference or a segment.'''

    return _KERNELS[edges.dtype].I_udv(edges, D_u, n_u, D_v, n_v, k)


def I_dudv(edges, D_u: float_array, n_u: float_array, D_v: float_array, n_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E \nabla u \cdot\mathbf{n} \overline{\nabla v \cdot \mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is either an arc of circunference or a segment.'''

    return _KERNELS[edges.dtype].I_dudv(edges, D_u, n_u, D_v, n_v, k)


def I_pw_v(edges, plane_wave: PlaneWave, D_v: float_array, n_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
.. math ::
        \int_E u_{\mathrm{inc}} \overline{v}\,\mathrm{d}\ell

    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is either an arc of circunference or a segment.'''

    return _KERNELS[edges.dtype].I_pw_v(edges, plane_wave, D_v, n_v, k)


def I_pw_dv(edges, plane_wave: PlaneWave, D_v: float_array, n_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is either an arc of circunference or a segment.'''

    return _KERNELS[edges.dtype].I_pw_dv(edges, plane_wave, D_v, n_v, k)


def I_dpw_v(edges, plane_wave: PlaneWave, D_v: float_array, n_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E \nabla u_{\mathrm{inc}}\cdot\mathbf{n} \overline{v}\,\mathrm{d}\ell

    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is either an arc of circunference or a segment.'''

    return _KERNELS[edges.dtype].I_dpw_v(edges, plane_wave, D_v, n_v, k)