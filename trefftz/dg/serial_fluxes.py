from numpy import dot
from trefftz.numpy_types import float_array
from enum import Enum
from typing import Protocol

from trefftz.dg.serial_kernels import I_uv, I_uv_arc, I_duv_arc, I_v_arc, I_dv_arc


JAC_ANGER_MODES = 80

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
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float, sign: SIGN) -> complex:
        ...


class SerialLocalKernel(Protocol):
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        ...

    def RHS(self, edge, d_psi: float_array, k: float) -> complex:
        ...


class SerialNonLocalKernel(Protocol):
    def LHS(self, edge_u, edge_v, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        ...


class NeumannFlux:
    '''Serial Neumann kernel'''
    def __init__(self, d_1: float):
        self.d_1 = d_1
    
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        d_1 = self.d_1
        d_m = d_psi
        d_n = d_phi

        N = edge["N"]

        return -1j*k*(1 + d_1 * dot(d_n, N))*dot(d_m, N)*I_uv(edge, d_phi, d_psi, k)


    def RHS(self, edge, d_psi: float_array, k: float) -> complex:
        raise NotImplementedError("Not implemented yet")



class UltraWeakFlux:
    '''Transmission flux for the UWVF'''
    def __init__(self, a: float, b: float):
        self.a = a 
        self.b = b
    
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float, sign: SIGN) -> complex:
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

        M = edge["M"]
        l = edge["l"]
        N = edge["N"]
        T = edge["T"]

        # I = -1j*k*l*((1/2+b*dot(d_n, N))*dot(d_m, N) + 1/2*dot(d_n, N) + a)*exp(1j*k*dot(d_n-d_m, M))*sinc(k*l/(2*pi)*dot(d_n-d_m, T))
        I = -1j*k*((1/2+b*dot(d_n, N))*dot(d_m, N) + 1/2*dot(d_n, N) + a)*I_uv(edge, d_n, d_m, k)

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


class DirichletFlux:
    '''Serial Dirichlet kernel'''
    def __init__(self, a: float, data = None):
        self.a = a
        self.data = data
    
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        a = self.a
        d_n = d_phi
        N = edge["N"]

        return -1j*k*(dot(d_n, N) + a)*I_uv(edge, d_phi, d_psi, k)
        
    def RHS(self, edge, d_psi: float_array, k: float) -> complex:
        raise NotImplementedError("Not implemented yet")


class CircularDirichletFlux:
    '''Serial Dirichlet kernel'''
    def __init__(self, a: float, data = None):
        self.a = a
        self.data = data
    
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        a = self.a

        return -I_duv_arc(edge, d_phi, d_psi, k) - 1j*k*a*I_uv_arc(edge, d_phi, d_psi, k)
        
    def RHS(self, edge, d_psi: float_array, k: float) -> complex:
        d_inc = self.data["d_inc"]
        a = self.a
        return I_dv_arc(edge, d_inc, d_psi, k) + 1j*a*k*I_v_arc(edge, d_inc, d_psi, k)


