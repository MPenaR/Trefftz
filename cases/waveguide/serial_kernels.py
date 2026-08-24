from trefftz.numpy_types import float_array
from numpy import conj, exp
from cases.waveguide.NeumannToDirichlet import Fdudn, Fu, ntd, Mode, ntd_mode, F_mode
from trefftz.dg.kernels.serial.linear_kernels import I_pw_v, I_pw_dv, PlaneWave


def I_Nuv(segment_u, segment_v, d_u: float_array, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:
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

# I could and SHOULD refactor it so I use the fact that u_inc = +- N(nabla(u_inc)) depending on the boundary you are. 
#any travelling mode is composed of two plane waves
def I_mode_v(segment_v, mode: Mode, d_v: float_array, k: float) -> complex:
    I = I_pw_v(segment_v, PlaneWave(mode.d_1, 0.5), d_v, k) + I_pw_v(segment_v, PlaneWave(mode.d_2, 0.5), d_v, k)
    return I


def I_mode_dv(segment_v, mode: Mode, d_v: float_array, k: float) -> complex:
    I = I_pw_dv(segment_v, PlaneWave(mode.d_1, 0.5), d_v, k) + I_pw_dv(segment_v, PlaneWave(mode.d_2, 0.5), d_v, k)
    return I


def I_Nmodedv(segment_v, mode: Mode, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:
    I = sum(ntd_mode(segment_v, mode, t)*conj(Fdudn(segment_v, d_v, k, H, t)) for t in range(0, NtD_modes))
    return I


def I_NmodeNv(segment_v, mode: Mode, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:
    I = sum(ntd_mode(segment_v, mode, t)*conj(ntd(segment_v, d_v, k, H, t)) for t in range(0, NtD_modes))
    return I


def I_Nmodev(segment_v, mode: Mode, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:
    I = sum(ntd_mode(segment_v, mode, t)*conj(Fu(segment_v, d_v, k, H, t)) for t in range(0, NtD_modes))
    return I


def I_modeNv(segment_v, mode: Mode, d_v: float_array, k: float, H: float,  NtD_modes: int) -> complex:
    I = sum(F_mode(segment_v, mode, t)*conj(ntd(segment_v, d_v, k, H, t)) for t in range(0, NtD_modes))
    return I
