from trefftz.numpy_types import float_array, complex_array
from numpy import outer, dot, exp, sinc, pi

def SoundHard(k: complex, l: float, N: float_array, T: float_array, M: float_array, d: float_array, d_d: float_array, d_1: float) -> complex_array:
    
    r"""
    Computes the block for an edge in a sound hard boundary.
    That is it computes the matrix :math:`\mathbf{M}=(M_{mn})` with:    
    
    .. math::
    
        M_{mn}=\boxed{(-ikl\left(1+d_{1}\mathbf{d}_{n}\cdot\mathbf{n}\right)\mathbf{d}_{m}\cdot\mathbf{n}e^{ik\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\cdot\mathbf{M}}\mathrm{sinc}\left(\frac{kl}{2\pi}\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\cdot\boldsymbol{\tau}\right)}


    Parameters:
    -----------

    - k : complex
        Wavenumber
    - l: float
        Length of the edge
    - N: float_array
        Normal of the edge
    - T: float_array
        Tangent of the edge
    - M: float_array
        Midpoint of the edge
    - d : float_array
        Set of directions
    - d_d : float_array
        Nd x Nd x 2 "Matrix" of differences of directions.
    - d_1 : float
        Stabilyzing parameter.
    """

    return -1j*k*l*outer(dot(d, N), (1 + d_1*dot(d, N)))*exp(1j*k*dot(d_d, M))*sinc(k*l/(2*pi)*dot(d_d, T))


def SoundSoft(k: complex, l: float, N: float_array, T: float_array, M: float_array, d: float_array, d_d: float_array, d_1: float) -> complex_array:
    pass


def Interior(k: complex, l: float, N: float_array, T: float_array, M: float_array, d: float_array, d_d: float_array, a: float, b: float) -> complex_array:
    pass
