from trefftz.numpy_types import float_array, complex_array
from numpy import conj
from cases.open_space.NeumannToDirichlet import Fu, Fdudn, ntd, ntd_inc, Fu_inc

def I_Nuv(arc_u, arc_v, D_u: float_array, D_v: float_array, k: float, NtD_modes: int, JA_modes: int) -> complex_array:
    R = arc_u["R"]
    I = R*sum(ntd(arc_u, D_u, k, t, JA_modes)[None, :]*conj(Fu(arc_v, D_v, k, t, JA_modes))[:, None] for t in range(-NtD_modes, NtD_modes+1))
    return I


def I_Nuincv(arc_v, d_inc: float_array, D_v: float_array, k: float, NtD_modes: int, JA_modes: int) -> complex_array:
    R = arc_v["R"]
    I = R*sum(ntd_inc(R, d_inc, k, t, JA_modes)*conj(Fu(arc_v, D_v, k, t, JA_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I



def I_uNv(arc_u, arc_v, D_u: float_array, D_v: float_array, k: float, NtD_modes: int, JA_modes: int) -> complex_array:
    R = arc_u["R"]
    I = R*sum(Fu(arc_u, D_u, k, t, JA_modes)[None, :]*conj(ntd(arc_v, D_v, k, t, JA_modes))[:, None] for t in range(-NtD_modes, NtD_modes+1))
    return I

def I_uincNv(arc_v, d_inc: float_array, D_v: float_array, k: float, NtD_modes: int, JA_modes: int) -> complex_array:
    R = arc_v["R"]
    I = R*sum(Fu_inc(R, d_inc, k, t, JA_modes)*conj(ntd(arc_v, D_v, k, t, JA_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I


def I_NuNv(arc_u, arc_v, D_u: float_array, D_v: float_array, k: float, NtD_modes: int, JA_modes: int) -> complex_array:
    R = arc_u["R"]
    I = R*sum(ntd(arc_u, D_u, k, t, JA_modes)[None, :]*conj(ntd(arc_v, D_v, k, t, JA_modes))[:, None] for t in range(-NtD_modes, NtD_modes+1))
    return I

def I_NuincNv(arc_v, d_inc: float_array, D_v: float_array, k: float, NtD_modes: int, JA_modes: int) -> complex_array:
    R = arc_v["R"]
    I = R*sum(ntd_inc(R, d_inc, k, t, JA_modes)*conj(ntd(arc_v, D_v, k, t, JA_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I

def I_Nudv(arc_u, arc_v, D_u: float_array, D_v: float_array, k: float, NtD_modes: int, JA_modes: int) -> complex_array:
    R = arc_u["R"]
    I = R*sum(ntd(arc_u, D_u, k, t, JA_modes)[None, :]*conj(Fdudn(arc_v, D_v, k, t, JA_modes))[:, None] for t in range(-NtD_modes, NtD_modes+1))
    return I


def I_Nuincdv(arc_v, d_inc: float_array, D_v: float_array, k: float, NtD_modes: int, JA_modes: int) -> complex_array:
    R = arc_v["R"]
    I = R*sum(ntd_inc(R, d_inc, k, t, JA_modes)*conj(Fdudn(arc_v, D_v, k, t, JA_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I