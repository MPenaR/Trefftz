from numpy import pi, exp, sinc, dot, array, cos
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
        self.direction = direction
        self.beta = beta(k, H, t)
        self.d_1 = array((direction*sqrt(1 - (t*pi/(k*H))**2), t*pi/(k*H)))
        self.d_2 = array((direction*sqrt(1 - (t*pi/(k*H))**2),-t*pi/(k*H)))


    def __call__(self, x: float, y: float):
        pm = self.direction
        beta = self.beta
        t = self.t
        H = self.H
        return exp(1j*pm*beta*x)*cos(t*pi*y/H)

def F_mode(segment, mode: Mode, t: int) -> complex:
    '''fourier coefficient of a mode restricted at region Sigma'''
    x = segment["P"][0]
    H = mode.H
    beta = mode.beta
    pm = mode.direction
    if mode.t != t:
        return 0
    else:
        if t == 0:
            return sqrt(H)*exp(pm*1j*beta*x)
        else:
            return sqrt(H/2)*exp(pm*1j*beta*x)

def F_dmodedn(segment, mode: Mode, t: int) -> complex:
    '''fourier coefficient of a mode restricted at region Sigma'''
    x = segment["P"][0]
    H = mode.H
    beta = mode.beta
    N = segment["N"]
    pm = mode.direction
    if mode.t != t:
        return 0
    else:
        if t == 0:
            return sqrt(H)*exp(pm*1j*beta*x)*pm*1j*beta*N[0]
        else:
            return sqrt(H/2)*exp(pm*1j*beta*x)*pm*1j*beta*N[0]



# #specific one
# def Fmode(segment, mode: Mode, t: int) -> complex:
#     H = mode.H
#     if mode.t != t:
#         return 0
#     else:
#         if t == 0:
#             return sqrt(H)
#         else:
#             return sqrt(H/2)



def ntd_mode(segment, mode: Mode, t: int) -> complex:
    beta = mode.beta
    return 1/(1j*beta)*F_dmodedn(segment, mode, t)
