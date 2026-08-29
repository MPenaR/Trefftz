from trefftz.numpy_types import float_array, complex_array
import numpy as np
from scipy.special import jn, jvp, hankel1, h1vp

def CircularDirichlet( X: float_array, Y: float_array, R: float, k: float, theta_inc: float, c: float_array = np.array([0., 0.]), M: int = 20, A: complex = 1. + 0j) -> complex_array:
    r = np.sqrt((X-c[0])**2 + (Y-c[1])**2)
    theta = np.atan2(Y - c[1], X - c[0])

    d = np.array([np.cos(theta_inc), np.sin(theta_inc)])

    r = np.expand_dims(r, 0)
    theta = np.expand_dims(theta, 0)
    p = np.expand_dims(np.arange(-M, M+1), list(np.arange(X.ndim, dtype=int)+1))
    u = - A*np.sum(1j**p*jn(p, k*R)*hankel1(p, k*r)/hankel1(p, k*R) * np.exp(1j*k*np.dot(d, c))*np.exp(1j*p*(theta-theta_inc)), 0)

    return np.where(r[0] > R, u, np.full_like(u, np.nan))


def nf_diel_cylinder_plane_wave( X: float_array, Y: float_array, k: float,
                          theta_inc: float, n: float, c: float_array, R: float,
                          M: int = 10, A: complex = 1+0j) -> complex_array:
    '''
    Scattered wave for a dielectric cylinder when irradiated with a plane wave
    Inputs:
    - X: x-coordinate of the points of evaluation
    - Y: y-coordinate of the points of evaluation
    - k: wavenumber
    - theta_inc: angle of propagation of the incident wave
      using the convention u_inc(x) = A exp(1j*k*d(theta_inc) dot x)
    - n: refraction index, square root of the relative electrical permittivity of the scatterer
    - c: location of the center of the cylinder
    - R: radius of the cylinder
    - M: number of modes uses in the series solution
    - A: Complex amplitude of the wave at (x,y)=(0,0) (default=1)
    returns
    - u: complex amplitude of the scattered field at every (x,y)
    '''

    r = np.sqrt((X-c[0])**2 + (Y-c[1])**2)
    theta = np.atan2(Y - c[1], X - c[0])
    r = np.expand_dims(r, 0)
    theta = np.expand_dims(theta, 0)

    d = np.array([np.cos(theta_inc), np.sin(theta_inc)])

    p = np.expand_dims(np.arange(-M, M+1), tuple(list( np.arange(X.ndim, dtype=int)+1)))  #set of "p" values in the "infinite" sum
    kR = k*R
    
    W = jn(p, n*kR)*h1vp(p, kR) - n*jvp(p, n*kR)*hankel1(p, kR)
    a = 1j**p * (n*jn(p, kR)*jvp(p, n*kR)-jn(p, n*kR)*jvp(p, kR)) / W
    b = 1j**p * (jn(p, kR)*h1vp(p, kR) - jvp(p, kR)*hankel1(p, kR)) / W
    u = A*np.exp(1j*k*np.dot(d, c))*np.sum(a*hankel1(p, k*r)*np.exp(1j*p*(theta-theta_inc)), 0)
    w = A*np.exp(1j*k*np.dot(d, c))*np.sum(b*     jn(p, n*k*r)*np.exp(1j*p*(theta-theta_inc)), 0)

    return np.where(r[0] > R, u, w)
