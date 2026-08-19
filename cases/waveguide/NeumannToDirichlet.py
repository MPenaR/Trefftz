from numpy import pi, exp, sinc, dot, array
from numpy.lib.scimath import sqrt
from trefftz.numpy_types import float_array
from enum import IntEnum

class PropagationDirection(IntEnum): 
    LEFT = -1
    RIGHT = +1

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


class Mode:
    def __init__(self, k: float, H: float, t: int, direction: PropagationDirection):
        self.k = k
        self.H = H
        self.t = t
        self.beta = beta(k, H, t)
        self.d_1 = array((direction*sqrt(1 - (t*pi/(k*H))**2), t*pi/(k*H)))
        self.d_2 = array((direction*sqrt(1 - (t*pi/(k*H))**2),-t*pi/(k*H)))



