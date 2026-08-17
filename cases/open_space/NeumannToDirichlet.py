from trefftz.numpy_types import float_array
from numpy import exp, atan2, sqrt, pi, conj
from scipy.special import jv, jvp, hankel1, h1vp

def Fu(arc, d: float_array, k: float, t: int, N_modes: int) -> complex:
    r'''
    Fourier coefficients of u
    '''

    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    phi = atan2(d[1], d[0])

    def primitive(theta: float) -> complex:
        I = 1j**t*exp(-1j*phi*t)*(jv(t, k*R)*theta + sum( 1j**p/(1j*p)*(jv(p+t, k*R)*exp( 1j*p*(theta - phi))
                                                              - (-1)**t*jv(p-t, k*R)*exp(-1j*p*(theta - phi))) for p in range(1, N_modes)))
        return I

    I = (primitive(theta=theta_2) - primitive(theta=theta_1))/sqrt(2*pi)

    return I


def Fu_inc(R: float, d_inc: float_array, k: float, t: int, N_modes: int) -> complex:
    r'''
    Fourier coefficients of u_inc (plane wave).
    '''

    theta_1 = 0
    theta_2 = 2*pi

    phi = atan2(d_inc[1], d_inc[0])

    def primitive(theta: float) -> complex:
        I = 1j**t*exp(-1j*phi*t)*(jv(t, k*R)*theta + sum( 1j**p/(1j*p)*(jv(p+t, k*R)*exp( 1j*p*(theta - phi))
                                                              - (-1)**t*jv(p-t, k*R)*exp(-1j*p*(theta - phi))) for p in range(1, N_modes)))
        return I

    I = (primitive(theta=theta_2) - primitive(theta=theta_1))/sqrt(2*pi)

    return I


# def Fu_inc2(arc, d_inc: float_array, k: float, t: int, N_modes: int) -> complex:
#     r'''
#     Fourier coefficients of u_inc (plane wave).
#     '''

#     theta_1 = 0
#     theta_2 = 2*pi
#     R = arc["R"]

#     phi = atan2(d_inc[1], d_inc[0])
#     I = 1j**t*exp(-1j*phi*t)*(jv(t, k*R)*2*pi + sum( 1j**p/(1j*p)*(jv(p+t, k*R)*exp( 1j*p*(theta - phi))
#                                                          - (-1)**t*jv(p-t, k*R)*exp(-1j*p*(theta - phi))) for p in range(1, N_modes)))
                              

#     I = (primitive(theta=theta_2) - primitive(theta=theta_1))/sqrt(2*pi)

#     return I



def Fdudn(edge, d: float_array, k: float, t: int, N_modes: int) -> complex:
    r'''
    Fourier coefficients of du/dn
    '''
    theta_1 = edge["theta_1"]
    theta_2 = edge["theta_2"]
    R = edge["R"]

    phi = atan2(d[1], d[0])

    def primitive(theta: float) -> complex:
        I = -1j*exp(-1j*t*phi)*1j**t*(jvp(t, k*R)*theta +
                                          sum( 1j**p/(1j*p)*(jvp(p+t, k*R)*exp( 1j*p*(theta - phi))
                                                    -(-1)**t*jvp(p-t, k*R)*exp(-1j*p*(theta - phi)) ) for p in range(1, N_modes)))
        return I

    I = 1j*k*(primitive(theta=theta_2) - primitive(theta=theta_1))/sqrt(2*pi)
    return I

def Fdu_incdn(R: float, d_inc: float_array, k: float, t: int, N_modes: int) -> complex:
    r'''
    Fourier coefficients of du/dn
    '''
    theta_1 = 0. 
    theta_2 = 2*pi

    phi = atan2(d_inc[1], d_inc[0])

    def primitive(theta: float) -> complex:
        I = -1j*exp(-1j*t*phi)*1j**t*(jvp(t, k*R)*theta +
                                          sum( 1j**p/(1j*p)*(jvp(p+t, k*R)*exp( 1j*p*(theta - phi))
                                                    -(-1)**t*jvp(p-t, k*R)*exp(-1j*p*(theta - phi)) ) for p in range(1, N_modes)))
        return I

    I = 1j*k*(primitive(theta=theta_2) - primitive(theta=theta_1))/sqrt(2*pi)
    return I



def NtD_coeff(edge, d: float_array, k: float, t: int, N_modes: int) -> complex:
    R = edge["R"]
    return 1/k*Fdudn(edge, d, k, t, N_modes)*hankel1(t, k*R)/h1vp(t, k*R)


def NtD_coeff_inc(R, d_inc: float_array, k: float, t: int, N_modes: int) -> complex:
    return 1/k*Fdu_incdn(R, d_inc, k, t, N_modes)*hankel1(t, k*R)/h1vp(t, k*R)


def I_Nuv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, NtD_modes: int, N_modes: int) -> complex:
    R = edge_u["R"]
    I = R*sum(NtD_coeff(edge_u, d_u, k, t, N_modes)*conj(Fu(edge_v, d_v, k, t, N_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I


def I_Nuincv(edge_v, d_inc: float_array, d_v: float_array, k: float, NtD_modes: int, N_modes: int) -> complex:
    R = edge_v["R"]
    I = R*sum(NtD_coeff_inc(R, d_inc, k, t, N_modes)*conj(Fu(edge_v, d_v, k, t, N_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I



def I_uNv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, NtD_modes: int, N_modes: int) -> complex:
    R = edge_u["R"]
    I = R*sum(Fu(edge_u, d_u, k, t, N_modes)*conj(NtD_coeff(edge_v, d_v, k, t, N_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I

def I_uincNv(edge_v, d_inc: float_array, d_v: float_array, k: float, NtD_modes: int, N_modes: int) -> complex:
    R = edge_v["R"]
    I = R*sum(Fu_inc(R, d_inc, k, t, N_modes)*conj(NtD_coeff(edge_v, d_v, k, t, N_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I


def I_NuNv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, NtD_modes: int, N_modes: int) -> complex:
    R = edge_u["R"]
    I = R*sum(NtD_coeff(edge_u, d_u, k, t, N_modes)*conj(NtD_coeff(edge_v, d_v, k, t, N_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I

def I_NuincNv(edge_v, d_inc: float_array, d_v: float_array, k: float, NtD_modes: int, N_modes: int) -> complex:
    R = edge_v["R"]
    I = R*sum(NtD_coeff_inc(R, d_inc, k, t, N_modes)*conj(NtD_coeff(edge_v, d_v, k, t, N_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I

def I_Nudv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, NtD_modes: int, N_modes: int) -> complex:
    R = edge_u["R"]
    I = R*sum(NtD_coeff(edge_u, d_u, k, t, N_modes)*conj(Fdudn(edge_v, d_v, k, t, N_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I


def I_Nuincdv(edge_v, d_inc: float_array, d_v: float_array, k: float, NtD_modes: int, N_modes: int) -> complex:
    R = edge_v["R"]
    I = R*sum(NtD_coeff_inc(R, d_inc, k, t, N_modes)*conj(Fdudn(edge_v, d_v, k, t, N_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I