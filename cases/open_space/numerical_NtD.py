import numpy as np 
from numpy import sqrt, exp, dot, trapezoid
from trefftz.numpy_types import float_array, complex_array
from scipy.special import hankel1, h1vp

def ntd(theta: float_array, du_dn: complex_array, k: float, R: float, t: int) -> complex:
    return 1/k*hankel1(t, k*R)/h1vp(t, k*R)*trapezoid(du_dn*exp(-1j*t*theta)/np.sqrt(2*np.pi), theta)

def NtD(theta: float_array, du_dn: complex_array, k: float, R: float, M: int) -> complex_array:
    return sum(ntd(theta, du_dn, k, R, t)*exp(1j*t*theta) for t in range(-M, M+1))


# def FourierCoefficient(theta: float_array, f: float_array, t: int) -> complex:
#     r'''Computes the t Fourier coefficient defined as:
    
#     .. math ::
#         \frac{1}{2\pi}*\int_0^{2\pi}f(\theta)e^{-it\theta}\,\mathrm{d}\theta
#     '''
#     return 1/sqrt(2*np.pi)*trapezoid(f*exp(-1j*t*theta), theta)


# def num_Fdudn(edge, d: float_array, k: float, R: float, t: int) -> complex:
#     r'''Computes the t Fourier coefficient of 
#     .. math ::
#     \nabla u\cdot \mathbf{n} 
#     supported at a single edge.
#     '''
#     P = edge["P"]
#     Q = edge["Q"]
#     theta_1 = np.atan2(P[1], P[0])
#     theta_2 = np.atan2(Q[1], Q[0])

#     theta = np.linspace(theta_1, theta_2, N_POINTS)
#     u_r = np.column_stack([np.cos(theta), np.sin(theta)])
#     N = u_r
#     x = R*u_r

#     u = exp(1j*k*dot(x, d))
#     du_dn = 1j*k*dot(N, d)*u

#     I = FourierCoefficient(theta, du_dn, t)

#     return I

# def num_Fu(edge, d: float_array, k: float, R: float, t: int) -> complex:
#     r'''Computes the t Fourier coefficient of 
#     .. math ::
#     \nabla u\cdot \mathbf{n} 
#     supported at a single edge.
#     '''
#     P = edge["P"]
#     Q = edge["Q"]
#     theta_1 = np.atan2(P[1], P[0])
#     theta_2 = np.atan2(Q[1], Q[0])

#     theta = np.linspace(theta_1, theta_2, N_POINTS)
#     u_r = np.column_stack([np.cos(theta), np.sin(theta)])
#     x = R*u_r

#     u = exp(1j*k*dot(x, d))
#     I = FourierCoefficient(theta, u, t)

#     return I

