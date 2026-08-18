from trefftz.numpy_types import float_array
from numpy import conj
from cases.waveguide.NeumanToDirichlet import Fdudn, Fu, ntd


def I_Nuv(segment_u, segment_v, d_u: float_array, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:
    # I = ntd(segment_u, d_u, k, H, 0)*conj(Fu(segment_v, d_v, k, H, 0))
    # I += sum(ntd(segment_u, d_u, k, H, 0)*conj(Fu(segment_v, d_v, k, H, 0)) for t in range(1, NtD_modes))
    I = sum(ntd(segment_u, d_u, k, H, t)*conj(Fu(segment_v, d_v, k, H, t)) for t in range(0, NtD_modes))

    return I

def I_uNv(segment_u, segment_v, d_u: float_array, d_v: float_array, k: float, H: float,  NtD_modes: int) -> complex:
    I = sum(Fu(segment_u, d_u, k, H, t)*conj(ntd(segment_v, d_v, k, H, t)) for t in range(0, NtD_modes))
    return I

def I_NuNv(segment_u, segment_v, d_u: float_array, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:
    I = sum(ntd(segment_u, d_u, k, H, t)*conj(ntd(segment_v, d_v, k, H, t)) for t in range(0, NtD_modes))
    return I

def I_Nudv(segment_u, segment_v, d_u: float_array, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:
    I = sum(ntd(segment_u, d_u, k, H, t)*conj(Fdudn(segment_v, d_v, k, H, t)) for t in range(0, NtD_modes))
    return I
