from trefftz.numpy_types import float_array
from numpy import conj
from cases.waveguide.NeumanToDirichlet import Fdudn, Fu, ntd


def I_Nuv(segment_u, segment_v, D_u: float_array, D_v: float_array, k: float, H: float, NtD_modes: int) -> complex:
    I = sum(ntd(segment_u, D_u, k, H, t)[None, :]*conj(Fu(segment_v, D_v, k, H, t))[:, None] for t in range(0, NtD_modes))
    return I

def I_uNv(segment_u, segment_v, D_u: float_array, D_v: float_array, k: float, H: float,  NtD_modes: int) -> complex:
    I = sum(Fu(segment_u, D_u, k, H, t)[None, :]*conj(ntd(segment_v, D_v, k, H, t))[:, None] for t in range(0, NtD_modes))
    return I

def I_NuNv(segment_u, segment_v, D_u: float_array, D_v: float_array, k: float, H: float, NtD_modes: int) -> complex:
    I = sum(ntd(segment_u, D_u, k, H, t)[None, :]*conj(ntd(segment_v, D_v, k, H, t))[:, None] for t in range(0, NtD_modes))
    return I

def I_Nudv(segment_u, segment_v, D_u: float_array, D_v: float_array, k: float, H: float, NtD_modes: int) -> complex:
    I = sum(ntd(segment_u, D_u, k, H, t)[None, :]*conj(Fdudn(segment_v, D_v, k, H, t))[:, None] for t in range(0, NtD_modes))
    return I
