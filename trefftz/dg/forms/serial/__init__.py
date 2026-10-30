from trefftz.numpy_types import float_array
from trefftz.mesh.core import edge_dtype, arc_dtype
from . import circular, linear
from trefftz.dg.exact import PlaneWave

_FORMS_LIBRARY = {
    edge_dtype: linear,
    arc_dtype: circular
}


def I_uv(edge, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is either an arc of circunference or a segment.'''

    return _FORMS_LIBRARY[edge.dtype].I_uv(edge, d_u, d_v, k)



def I_duv(edge, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E \nabla u\cdot\mathbf{n} \overline{v}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is either an arc of circunference or a segment.'''

    return _FORMS_LIBRARY[edge.dtype].I_duv(edge, d_u, d_v, k)


def I_udv(edge, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u \overline{\nabla v \cdot \mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is either an arc of circunference or a segment.'''

    return _FORMS_LIBRARY[edge.dtype].I_udv(edge, d_u, d_v, k)


def I_dudv(edge, d_u: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E \nabla u \cdot\mathbf{n} \overline{\nabla v \cdot \mathbf{n}}\,\mathrm{d}\ell

    where $u$ and $v$ are plane waves and $E$ is either an arc of circunference or a segment.'''

    return _FORMS_LIBRARY[edge.dtype].I_dudv(edge, d_u, d_v, k)


def I_uincv(edge, d_inc: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
.. math ::
        \int_E u_{\mathrm{inc}} \overline{v}\,\mathrm{d}\ell

    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is either an arc of circunference or a segment.'''

    return _FORMS_LIBRARY[edge.dtype].I_uincv(edge, d_inc, d_v, k)


def I_uincdv(edge, d_inc: float_array, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is either an arc of circunference or a segment.'''

    return _FORMS_LIBRARY[edge.dtype].I_uincdv(edge, d_inc, d_v, k)


def I_pw_v(edge, plane_wave: PlaneWave, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
.. math ::
        \int_E u_{\mathrm{inc}} \overline{v}\,\mathrm{d}\ell

    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is either an arc of circunference or a segment.'''

    return _FORMS_LIBRARY[edge.dtype].I_pw_v(edge, plane_wave, d_v, k)


def I_pw_dv(edge, plane_wave: PlaneWave, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E u_{\mathrm{inc}} \overline{\nabla v\cdot\mathbf{n}}\,\mathrm{d}\ell

    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is either an arc of circunference or a segment.'''

    return _FORMS_LIBRARY[edge.dtype].I_pw_dv(edge, plane_wave, d_v, k)


def I_dpw_v(edge, plane_wave: PlaneWave, d_v: float_array, k: float) -> complex:
    r'''Computes the integral:
    .. math ::
        \int_E \nabla u_{\mathrm{inc}}\cdot\mathbf{n} \overline{v}\,\mathrm{d}\ell

    where $u_inc$ is an incident plane wave and $v$ is a plane wave and $E$ is either an arc of circunference or a segment.'''

    return _FORMS_LIBRARY[edge.dtype].I_dpw_v(edge, plane_wave, d_v, k)