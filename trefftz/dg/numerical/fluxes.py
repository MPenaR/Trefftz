from numpy import dot, conj, linspace, outer, exp, column_stack, cos, sin
from trefftz.numpy_types import float_array
from trefftz.dg.serial_fluxes import SIGN
from numpy import trapezoid as Int

N_POINTS = 10**5

class NeumannFlux:
    '''Serial Neumann kernel'''
    def __init__(self, d_1: float):
        self.d_1 = d_1
    
    def LHS(self, edge, d_u: float_array, d_v: float_array, k: float) -> complex:
        d_1 = self.d_1

        N = edge["N"]
        l = edge["l"]
        P = edge["P"]
        T = edge["T"]
        
        t = linspace(0., 1., N_POINTS)
        x = P + outer(t, T)*l
        u = exp(1j*k*dot(x, d_u))
        du_dn = 1j*k*dot(N, d_u)*u
        v = exp(1j*k*dot(x, d_v))
        dv_dn = 1j*k*dot(N, d_v)*v

        return Int((u + d_1/(1j*k)*du_dn)*conj(dv_dn), t)*l


    def RHS(self, edge, d_v: float_array, k: float) -> complex:
        raise NotImplementedError("Not implemented yet")

def avg(w):
    return 0.5*w

def jump(w, sign: SIGN):
    return (-1)**sign.value[0]*w


class UltraWeakFlux:
    '''Transmission flux for the UWVF'''
    def __init__(self, a: float, b: float):
        self.a = a 
        self.b = b
    
    def LHS(self, edge, d_u: float_array, d_v: float_array, k: float, sign: SIGN) -> complex:
        a = self.a
        b = self.b

        N = edge["N"]
        l = edge["l"]
        P = edge["P"]
        T = edge["T"]
        
        t = linspace(0., 1., N_POINTS)
        x = P + outer(t, T)*l
        u = exp(1j*k*dot(x, d_u))
        du_dn = 1j*k*dot(N, d_u)*u
        v = exp(1j*k*dot(x, d_v))
        dv_dn = 1j*k*dot(N, d_v)*v

        I = Int((avg(u) + b/(1j*k)*jump(du_dn, sign))*conj(jump(dv_dn, sign)) - (a*1j*k*jump(u, sign) + avg(du_dn))*conj(jump(v, sign)), t)*l
        return I


class DirichletFlux:
    '''Serial Dirichlet kernel'''
    def __init__(self, a: float, data = None):
        self.a = a
        self.data = data
    
    def LHS(self, edge, d_u: float_array, d_v: float_array, k: float) -> complex:
        a = self.a

        N = edge["N"]
        l = edge["l"]
        P = edge["P"]
        T = edge["T"]
        
        t = linspace(0., 1., N_POINTS)
        x = P + outer(t, T)*l
        u = exp(1j*k*dot(x, d_u))
        du_dn = 1j*k*dot(N, d_u)*u
        v = exp(1j*k*dot(x, d_v))

        return -Int((du_dn + 1j*a*k*u)*conj(v), t)*l
        
    def RHS(self, edge, d_v: float_array, k: float) -> complex:
        d_inc = self.data["d_inc"]
        a = self.a

        N = edge["N"]
        l = edge["l"]
        P = edge["P"]
        T = edge["T"]
        
        t = linspace(0., 1., N_POINTS)
        x = P + outer(t, T)*l

        v = exp(1j*k*dot(x, d_v))
        dv_dn = 1j*k*dot(N, d_v)*v

        u_inc = -exp(1j*k*dot(x, d_inc))

        return -Int(u_inc*conj(dv_dn - 1j*a*k*v), t)*l

class CircularDirichletFlux:
    '''Serial Dirichlet kernel'''
    def __init__(self, a: float, data = None):
        self.a = a
        self.data = data
    
    def LHS(self, edge, d_u: float_array, d_v: float_array, k: float) -> complex:
        a = self.a

        R = edge["R"]

        theta_1 = edge["theta_1"]
        theta_2 = edge["theta_2"]
        theta = linspace(theta_1, theta_2, N_POINTS)

        u_r = column_stack([cos(theta), sin(theta)])
        x = R*u_r
        N = u_r
        
        u = exp(1j*k*dot(x, d_u))
        du_dn = 1j*k*dot(N, d_u)*u
        v = exp(1j*k*dot(x, d_v))

        return -Int((du_dn + 1j*a*k*u)*conj(v), theta)*R
        
    def RHS(self, edge, d_v: float_array, k: float) -> complex:
        d_inc = self.data["d_inc"]
        a = self.a

        R = edge["R"]

        theta_1 = edge["theta_1"]
        theta_2 = edge["theta_2"]
        theta = linspace(theta_1, theta_2, N_POINTS)

        u_r = column_stack([cos(theta), sin(theta)])
        x = R*u_r
        N = u_r
        
        v = exp(1j*k*dot(x, d_v))
        dv_dn = 1j*k*dot(N, d_v)*v

        u_inc = -exp(1j*k*dot(x, d_inc))

        return -Int(u_inc*conj(dv_dn - 1j*a*k*v), theta)*R
