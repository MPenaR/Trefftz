from numpy import pi, exp, sinc, dot
from numpy.lib.scimath import sqrt
from trefftz.numpy_types import float_array

def Fu(edge, d: float_array, k: float, H: float, t: int) -> complex:
    r'''
    "t" cosine coefficient of u
    '''

    M = edge["M"]
    l = edge["l"]
    if t == 0:
        I =  l/sqrt(H)*exp(1j*k*dot(d, M))*sinc(k*l/(2*pi)*d[...,1])
    else: 
        I = l/sqrt(2*H)*exp(1j*k*dot(d, M))*(exp( 1j*t*pi*M[1]/H)*sinc(l/(2*H)*(t + k*H*d[..., 1]/pi)) +
                                             exp(-1j*t*pi*M[1]/H)*sinc(l/(2*H)*(t - k*H*d[..., 1]/pi)))
    return I


def Fdudn(edge, d: float_array, k: float, H: float,  t: int) -> complex:
    r'''
    "t" cosine coefficient of du/dn
    '''
    N = edge["N"]
    return 1j*k*dot(d, N)*Fu(edge, d, k, H, t)


def beta(k: float, H: float, t: int) -> complex | float:
    return sqrt(k**2 - (t*pi/H)**2)


def ntd(edge, d: float_array, k: float, H: float, t: int) -> complex:
    return 1/(1j*beta(k, H, t))*Fdudn(edge, d, k, H, t)