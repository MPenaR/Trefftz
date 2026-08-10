from numpy import pi, exp, conj, sinc, dot
from numpy.lib.scimath import sqrt
from trefftz.numpy_types import float_array

def Fu(edge, d: float_array, k: float, H: float, t: int) -> complex:
    r'''
    "t" cosine coefficient of u
    '''

    M = edge["M"]
    l = edge["l"]
    if t == 0:
        I =  l/sqrt(H)*exp(1j*k*dot(d, M))*sinc(k*l/(2*pi)*d[1])
    else: 
        I = l/sqrt(2*H)*exp(1j*k*dot(d, M))*(exp( 1j*t*pi*M[1]/H)*sinc(l/(2*H)*(t + k*H*d[1]/pi)) +
                                             exp(-1j*t*pi*M[1]/H)*sinc(l/(2*H)*(t - k*H*d[1]/pi)))
    return I


def Fdudn(edge, d: float_array, k: float, H: float,  t: int) -> complex:
    r'''
    "t" cosine coefficient of du/dn
    '''
    N = edge["N"]
    return 1j*k*dot(d, N)*Fu(edge, d, k, H, t)


def beta(k: float, H: float, t: int) -> complex | float:
    return sqrt(k**2 - (t*pi/H)**2)

def NtD_coeff(edge, d: float_array, k: float, H: float, t: int) -> complex:
    return 1/(1j*beta(k, H, t))*Fdudn(edge, d, k, H, t)


def I_Nuv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:
    # I = NtD_coeff(edge_u, d_u, k, H, 0)*conj(Fu(edge_v, d_v, k, H, 0))
    # I += sum(NtD_coeff(edge_u, d_u, k, H, 0)*conj(Fu(edge_v, d_v, k, H, 0)) for t in range(1, NtD_modes))
    I = sum(NtD_coeff(edge_u, d_u, k, H, t)*conj(Fu(edge_v, d_v, k, H, t)) for t in range(0, NtD_modes))

    return I

def I_uNv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, H: float,  NtD_modes: int) -> complex:
    I = sum(Fu(edge_u, d_u, k, H, t)*conj(NtD_coeff(edge_v, d_v, k, H, t)) for t in range(0, NtD_modes))
    return I

def I_NuNv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:
    I = sum(NtD_coeff(edge_u, d_u, k, H, t)*conj(NtD_coeff(edge_v, d_v, k, H, t)) for t in range(0, NtD_modes))
    return I

def I_Nudv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:
    I = sum(NtD_coeff(edge_u, d_u, k, H, t)*conj(Fdudn(edge_v, d_v, k, H, t)) for t in range(0, NtD_modes))
    return I
