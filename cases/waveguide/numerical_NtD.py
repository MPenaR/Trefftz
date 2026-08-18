import numpy as np
from numpy import pi, cos, sqrt
from numpy import trapezoid
from trefftz.numpy_types import complex_array, float_array
from cases.waveguide.NeumannToDirichlet import beta


def ntd(y: float_array, du_dn: complex_array, k: float, H: float, t: int) -> complex:
    if t == 0:
        fdudn =  trapezoid(du_dn*1/np.sqrt(H), y )
    else:
        fdudn = trapezoid(du_dn*cos(t*pi*y/H)/np.sqrt(H/2), y )
    return 1/(1j*beta(k, H, t))*fdudn


def NtD(y: float_array, du_dn: complex_array, k: float, H: float, NtD_modes: int):

    u = ntd(y, du_dn, k, H, 0)/sqrt(H)*np.ones_like(y) + sum([ ntd(y, du_dn, k, H, t)*cos(t*pi*y/H)/np.sqrt(H/2) for t in range(1, NtD_modes)])
    return u
