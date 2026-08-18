from trefftz.numpy_types import float_array
from numpy import exp, atan2, sqrt, pi, conj
from scipy.special import jv, jvp, hankel1, h1vp

def Fu(arc, d: float_array, k: float, t: int, JA_modes: int) -> complex:
    r'''
    Fourier coefficients of u
    '''

    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    phi = atan2(d[..., 1], d[..., 0])

    def primitive(theta: float) -> complex:
        I = 1j**t*exp(-1j*phi*t)*(jv(t, k*R)*theta + sum( 1j**p/(1j*p)*(jv(p+t, k*R)*exp( 1j*p*(theta - phi))
                                                              - (-1)**t*jv(p-t, k*R)*exp(-1j*p*(theta - phi))) for p in range(1, JA_modes)))
        return I

    I = (primitive(theta=theta_2) - primitive(theta=theta_1))/sqrt(2*pi)

    return I


def Fu_inc(R: float, d_inc: float_array, k: float, t: int, JA_modes: int) -> complex:
    r'''
    Fourier coefficients of u_inc (plane wave).
    '''

    theta_1 = 0
    theta_2 = 2*pi

    phi = atan2(d_inc[1], d_inc[0])

    def primitive(theta: float) -> complex:
        I = 1j**t*exp(-1j*phi*t)*(jv(t, k*R)*theta + sum( 1j**p/(1j*p)*(jv(p+t, k*R)*exp( 1j*p*(theta - phi))
                                                              - (-1)**t*jv(p-t, k*R)*exp(-1j*p*(theta - phi))) for p in range(1, JA_modes)))
        return I

    I = (primitive(theta=theta_2) - primitive(theta=theta_1))/sqrt(2*pi)

    return I


def Fdudn(arc, d: float_array, k: float, t: int, JA_modes: int) -> complex:
    r'''
    Fourier coefficients of du/dn
    '''
    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    phi = atan2(d[..., 1], d[..., 0])

    def primitive(theta: float) -> complex:
        I = -1j*exp(-1j*t*phi)*1j**t*(jvp(t, k*R)*theta +
                                          sum( 1j**p/(1j*p)*(jvp(p+t, k*R)*exp( 1j*p*(theta - phi))
                                                    -(-1)**t*jvp(p-t, k*R)*exp(-1j*p*(theta - phi)) ) for p in range(1, JA_modes)))
        return I

    I = 1j*k*(primitive(theta=theta_2) - primitive(theta=theta_1))/sqrt(2*pi)
    return I

def Fdu_incdn(R: float, d_inc: float_array, k: float, t: int, JA_modes: int) -> complex:
    r'''
    Fourier coefficients of du/dn
    '''
    theta_1 = 0. 
    theta_2 = 2*pi

    phi = atan2(d_inc[1], d_inc[0])

    def primitive(theta: float) -> complex:
        I = -1j*exp(-1j*t*phi)*1j**t*(jvp(t, k*R)*theta +
                                          sum( 1j**p/(1j*p)*(jvp(p+t, k*R)*exp( 1j*p*(theta - phi))
                                                    -(-1)**t*jvp(p-t, k*R)*exp(-1j*p*(theta - phi)) ) for p in range(1, JA_modes)))
        return I

    I = 1j*k*(primitive(theta=theta_2) - primitive(theta=theta_1))/sqrt(2*pi)
    return I



def ntd(arc, d: float_array, k: float, t: int, JA_modes: int) -> complex:
    R = arc["R"]
    return 1/k*Fdudn(arc, d, k, t, JA_modes)*hankel1(t, k*R)/h1vp(t, k*R)


def ntd_inc(R, d_inc: float_array, k: float, t: int, JA_modes: int) -> complex:
    return 1/k*Fdu_incdn(R, d_inc, k, t, JA_modes)*hankel1(t, k*R)/h1vp(t, k*R)