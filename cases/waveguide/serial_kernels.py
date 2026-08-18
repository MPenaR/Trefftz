from trefftz.numpy_types import float_array
from numpy import conj
from cases.waveguide.NeumanToDirichlet import Fdudn, Fu, ntd


def I_Nuv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:
    # I = ntd(edge_u, d_u, k, H, 0)*conj(Fu(edge_v, d_v, k, H, 0))
    # I += sum(ntd(edge_u, d_u, k, H, 0)*conj(Fu(edge_v, d_v, k, H, 0)) for t in range(1, NtD_modes))
    I = sum(ntd(edge_u, d_u, k, H, t)*conj(Fu(edge_v, d_v, k, H, t)) for t in range(0, NtD_modes))

    return I

def I_uNv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, H: float,  NtD_modes: int) -> complex:
    I = sum(Fu(edge_u, d_u, k, H, t)*conj(ntd(edge_v, d_v, k, H, t)) for t in range(0, NtD_modes))
    return I

def I_NuNv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:
    I = sum(ntd(edge_u, d_u, k, H, t)*conj(ntd(edge_v, d_v, k, H, t)) for t in range(0, NtD_modes))
    return I

def I_Nudv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:
    I = sum(ntd(edge_u, d_u, k, H, t)*conj(Fdudn(edge_v, d_v, k, H, t)) for t in range(0, NtD_modes))
    return I
