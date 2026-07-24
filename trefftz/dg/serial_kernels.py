from numpy import dot, exp, sinc, pi, conj
from numpy.lib.scimath import sqrt
from typing import NamedTuple
from trefftz.numpy_types import float_array
from enum import Enum
from typing import Protocol


class Edge(NamedTuple):
    M: float_array
    l: float_array
    N: float_array
    T: float_array


class SIGN(Enum):
    '''Sign for the transmission kernel
    where PP (plus plus) stands for both trial
    and test function coming from the plus triangle
    whereas PM stands for the test function (psi) coming
    from the + triangle and the trial function (phi) from
    the - one.'''
    PP = (0, 0)
    PM = (0, 1)
    MP = (1, 0)
    MM = (1, 1)


class SerialTransmissionKernel(Protocol):
    def LHS(self, edge: Edge, d_phi: float_array, d_psi: float_array, k: float, sign: SIGN) -> complex:
        ...


class SerialLocalKernel(Protocol):
    def LHS(self, edge: Edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        ...

    def RHS(self, edge: Edge, d_psi: float_array, k: float) -> complex:
        ...


class SerialNonLocalKernel(Protocol):
    def LHS(self, edge_1: Edge, edge_2: Edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        ...


class NeumannKernel:
    '''Serial Neumann kernel'''
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


class DirichletKernel:
    '''Serial Dirichlet kernel'''
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
        raise NotImplementedError("Not implemented yet")
        
    def RHS(self, edge: Edge, d_psi: float_array, k: float) -> complex:
        raise NotImplementedError("Not implemented yet")



class UltraWeakKernel:
    '''Transmission kernel for the UWVF'''
    def __init__(self, a: float, b: float):
        self.a = a 
        self.b = b
    
    def LHS(self, edge: Edge, d_phi: float_array, d_psi: float_array, k: float, sign: SIGN) -> complex:
        match sign:
            case SIGN.PP:
                a = self.a
                b = self.b
            case SIGN.PM:
                a = self.a
                b = self.b
            case SIGN.MP:
                a = -self.a
                b = -self.b
            case SIGN.MM:
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
            case SIGN.PP:
                I = I
            case SIGN.PM:
                I = -I
            case SIGN.MP:
                I = I
            case SIGN.MM:
                I = -I
        return I