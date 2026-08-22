from trefftz.numpy_types import float_array, complex_array
from trefftz.mesh.core import edge_dtype, arc_dtype
from trefftz.dg.kernels.block import arc_kernels, linear_kernels

_KERNELS = {
    edge_dtype: linear_kernels,
    arc_dtype: arc_kernels
}


def I_uv(edge, D_u: float_array, D_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is either an arc of circunference or a segment.'''

    return _KERNELS[edge.dtype].I_uv(edge, D_u, D_v, k)



def I_duv(edge, D_u: float_array, D_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E \nabla u\cdot\mathbf{n} \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is either an arc of circunference or a segment.'''

    return _KERNELS[edge.dtype].I_duv(edge, D_u, D_v, k)


def I_udv(edge, D_u: float_array, D_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{\nabla v \cdot \mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is either an arc of circunference or a segment.'''

    return _KERNELS[edge.dtype].I_udv(edge, D_u, D_v, k)


def I_dudv(edge, D_u: float_array, D_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E \nabla u \cdot\mathbf{n} \overline{\nabla v \cdot \mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is either an arc of circunference or a segment.'''

    return _KERNELS[edge.dtype].I_dudv(edge, D_u, D_v, k)


def I_uincv(edge, d_inc: float_array, D_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
.. math ::
        \int_E u_{\mathrm{inc}} \overline{v}\,\mathrm{d}\ell

    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is either an arc of circunference or a segment.'''

    return _KERNELS[edge.dtype].I_uincv(edge, d_inc, D_v, k)


def I_uincdv(edge, d_inc: float_array, D_v: float_array, k: float) -> complex_array:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is either an arc of circunference or a segment.'''

    return _KERNELS[edge.dtype].I_uincdv(edge, d_inc, D_v, k)