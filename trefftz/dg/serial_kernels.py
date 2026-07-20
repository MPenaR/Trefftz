from numpy import dot, exp, sinc, pi, conj
from numpy.lib.scimath import sqrt
from typing import NamedTuple
from trefftz.numpy_types import float_array
from enum import Enum, auto

class Edge(NamedTuple):
    M: float_array
    l: float_array
    N: float_array
    T: float_array


class TT(Enum):
    '''sentinel for the transmission kernel
    where PP (plus plus) stands for both trial
    and test function coming from the plus triangle
    whereas PM stands for the test function (psi) coming
    from the + triangle and the trial function (phi) from
    the - one.'''
    PP = auto()
    PM = auto()
    MP = auto()
    MM = auto()

class SoundHardKernel:
    '''Serial SoundHard kernel'''
    def __init__(self, d_1: float):
        self.d_1 = d_1
    
    def LHS(self, edge: Edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        d_1 = self.d_1
        d_m = d_psi
        d_n = d_phi

        M = edge.M
        l = edge.l
        N = edge.N
        T = edge.T

        return -1j*k*l*(1 + d_1 * dot(d_n, N))*dot(d_m, N)*exp(1j*k*dot(d_n - d_m, M)) * sinc(k*l/(2*pi)*dot(d_n-d_m, T))
    
    def RHS(self, edge: Edge, d_psi: float_array, k: float) -> complex:
        raise NotImplementedError("Not implemented yet")


class UltraWeakKernel:
    '''Transmission kernel for the UWVF'''
    def __init__(self, a: float, b: float):
        self.a = a 
        self.b = b
    
    def LHS(self, edge: Edge, d_phi: float_array, d_psi: float_array, k: float, sign: TT) -> complex:
        match sign:
            case TT.PP:
                a = self.a
                b = self.b
            case TT.PM:
                a = -self.a
                b = -self.b
            case TT.MP:
                a = self.a
                b = self.b
            case TT.MM:
                a = -self.a
                b = -self.b

        d_m = d_psi
        d_n = d_phi

        k_n = k
        k_m = k

        M = edge.M
        l = edge.l
        N = edge.N
        T = edge.T

        # I = -1j*l/2*(2*a*k + k_n*dot(d_n, N) + k_m*dot(d_m, N) + 2*b/k*k_n*dot(d_n, N)*k_m*dot(d_m, N))*exp(1j*dot(k_n*d_n - k_m*d_m, M))*sinc(l/(2*pi)*dot(k_n*d_n - k_m*d_m,T))
        # I = -1j*k*l*( 1/2*( k_n/k*dot(d_n, N) + k_m/k*dot(d_m, N)) + a + b*k_n/k*k_m/k*dot(d_n, N)*dot(d_m, N))*exp(1j*dot(k_n*d_n - k_m*d_m, M))*sinc(l/(2*pi)*dot(k_n*d_n - k_m*d_m,T))
        I = -1j*k*l*((1/2+b*dot(d_n, N))*dot(d_m, N) + 1/2*dot(d_n, N) + a)*exp(1j*k*dot(d_n-d_m, M))*sinc(k*l/(2*pi)*dot(d_n-d_m, T))
        match sign:
            case TT.PP:
                I = I
            case TT.PM:
                I = I
            case TT.MP:
                I = -I
            case TT.MM:
                I = -I
        return I
class NtDLocal:
    def __init__(self, R: float, d_2: float, n: int, H: float):
        self.R = R
        self.mode_n = n
        self.H = H
        self.d_2 = d_2
    
    def LHS(self, edge: Edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        r"""
        Computes the flux on a radiating boundary with respect to the degrees
        of freedom from the same cell, that is:

        TODO: it is assuming that the radiating boundary consists of a vertical segment. This should be easy to generalize.
        
        .. math::
        
            -\int_{E}\left(d_{2}ik\Phi_n(\mathbf{x})+\nabla \Phi_n(\mathbf{x})\cdot\mathbf{n}\right)\overline{\Psi_n(\mathbf{x})}\,\mathrm{d}S_{\mathbf{x}}

        which can be computed as:

        .. math::
        
            \boxed{-ikl\left(d_{2}+\mathbf{d}_{n}\cdot\mathbf{n}\right)e^{ik\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\cdot\mathbf{M}}\mathrm{sinc}\left(\frac{kl}{2\pi}\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\mathbf{j}\right)}

        Parameters
        ----------
        phi : Function
            Trial function.
        psi : Function
            Test function.
        k : float
            Wave number.
        edge : Edge
            Edge parameters.
        d_2 : float
            Stabilyzing parameter.

        Returns
        -------
        I : complex
            The integral.
        """

        d_2 = self.d_2
        d_n = d_phi
        d_m = d_psi


        M = edge.M
        l = edge.l
        N = edge.N
        T = edge.T

        return -1j*k*l*(d_2 + dot(d_n, N))*exp(1j*k*dot(d_n - d_m, M))*sinc(k*l/(2*pi)*dot(d_n-d_m, T))

    def RHS(self, edge: Edge, d_psi: float_array, k: float) -> complex:

        t = self.mode_n
        d = d_psi
        d_2 = self.d_2
        H = self.H
        
        M = edge.M
        l = edge.l
        N = edge.N
        T = edge.T

        M_x, _ = M  # I dont like it, it still assumes horizontal waveguide

        beta = sqrt(k**2 - (t*pi/H)**2)

        if t == 0:
            I = 2*1j*k*l*(dot(d, N) - d_2*(1-dot(d, N)))*exp(1j*k*(M_x - dot(d, M)))*sinc(k*l/(2*pi)*dot(d, T))
        else:
            I = 1j*k*l*(dot(d, N) - d_2*(1-k/conj(beta)*dot(d, N)))*exp(1j*(beta*M_x - k*dot(d, M)))*( exp( 1j*t*pi*dot(M,T)/H)*sinc(k*l/(2*pi)*dot(d, T) - t*(l/(2*H)))+
                                                                                                       exp(-1j*t*pi*dot(M,T)/H)*sinc(k*l/(2*pi)*dot(d, T) + t*(l/(2*H))))
        return I




class WaveguideNtD_nonlocal:
    def __init__(self, R: float, d_2: float, M: int, H: float): 
        self.R = R
        self.d_2 = d_2 
        self.M = M 
        self.H = H

    def LHS(self, edge_1: Edge, edge_2:Edge, d_phi: float_array, d_psi: float_array, k : float) -> complex:
        r"""
        Computes the flux on a radiating boundary with respect to the degrees
        of freedom from another cell, that is:

        TODO: it is assuming that the radiating boundary consists of a vertical segment. This should be easy to generalize.
        
        Parameters
        ----------
        phi : Function
            Trial function.
        psi : Function
            Test function.
        k : float
            Wave number.
        edge_u : Edge
            Edge of the triangle associated to the trial function.
        edge_v : Edge
            Edge of the triangle associated to the test function.
        d_2 : float
            Stabilyzing parameter.
        NtD_modes : int
            Number of modes for the approximation of the NtD map.
        H : float
            height of the waveguide. 

        Returns
        -------
        I : complex
            The integral.

        
        """

        d_2 = self.d_2
        H = self.H
        M = self.M

        M_u = edge_1.M
        l_u = edge_1.l

        M_v = edge_2.M
        l_v = edge_2.l


        N = edge_1.N
        # T = edge_1.T

        d_n = d_phi
        d_m = d_psi

        I1 = -1j*k*H*dot(d_n,N)*dot(d_m,N)*d_2*exp(1j*k*(dot(d_n,M_u) - dot(d_m,M_v)))*l_u/H*l_v/H*(
            sinc(k*l_u/(2*pi)*d_n[1])*sinc(k*l_v/(2*pi)*d_m[1]) + 1/2*sum( [ k**2 / abs(sqrt(complex(k**2 - (s*pi/H)**2)))**2 * (
            exp( 1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] + s*l_u/(2*H)) + exp(-1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] - s*l_u/(2*H)) ) *(
            exp(-1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] + s*l_v/(2*H)) + exp( 1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] - s*l_v/(2*H)) )
            for s in range(1, M)]) )
        
        I2 = -1j*k*H*dot(d_n,N)*(dot(d_m,N)-d_2)*exp(1j*k*(dot(d_n,M_u) - dot(d_m,M_v)))*l_u/H*l_v/H*(
            sinc(k*l_u/(2*pi)*d_n[1])*sinc(k*l_v/(2*pi)*d_m[1]) + 1/2*sum( [ k / sqrt(complex(k**2 - (s*pi/H)**2)) * (
            exp( 1j*s*pi*M_u[1]/H)*sinc(k*l_u/(2*pi)*d_n[1] + s*l_u/(2*H)) + exp(-1j*s*pi*M_u[1]/H)*sinc(k*l_u/(2*pi)*d_n[1] - s*l_u/(2*H)) ) *(
            exp(-1j*s*pi*M_v[1]/H)*sinc(k*l_v/(2*pi)*d_m[1] + s*l_v/(2*H)) + exp( 1j*s*pi*M_v[1]/H)*sinc(k*l_v/(2*pi)*d_m[1] - s*l_v/(2*H)) )
            for s in range(1, M)]) )
        
        I3 = 1j*k*H*dot(d_m,N)*d_2*exp(1j*k*(dot(d_n,M_u) - dot(d_m,M_v)))*l_u/H*l_v/H*(
            sinc(k*l_u/(2*pi)*d_n[1])*sinc(k*l_v/(2*pi)*d_m[1]) + 1/2*sum( [ k / conj(sqrt(complex(k**2 - (s*pi/H)**2))) * (
            exp( 1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] + s*l_u/(2*H)) + exp(-1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] - s*l_u/(2*H)) ) *(
            exp(-1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] + s*l_v/(2*H)) + exp( 1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] - s*l_v/(2*H)) )
            for s in range(1, M)]) )

        return  I1 + I2 + I3
