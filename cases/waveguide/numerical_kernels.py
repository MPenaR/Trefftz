from numpy import linspace, outer, sin, cos, pi, exp, dot, conj, isclose, array
from numpy.lib.scimath import sqrt
from numpy.linalg import norm
from numpy import trapezoid as Int
import numpy as np
from scipy.special import hankel1, h1vp
from trefftz.numpy_types import float_array, complex_array

N_POINTS = int(1E5)
NtD_MODES = 20



class numerical_NtDLocal:
    def __init__(self, R: float, d_2: float, n: int):
        self.R = R
        self.mode_n = n
        self.d_2 = d_2

    
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        d2 = self.d_2
        d_n = d_phi
        d_m = d_psi

        P = edge["P"]

        t = np.linspace(0, 1, N_POINTS)
        N = edge["N"]
        T = edge["T"]
        l = edge["l"]
        x = P + np.outer(t,T)*l
  
        u = exp(1j*k*dot(x, d_n))
        v = exp(1j*k*dot(x, d_m))
        du_dn = 1j*k*dot(N, d_n)*u
     
        I_uv = Int(u*conj(v), t)*l  # int -i*k*d_2*u*conj(v) dl 
        I_duv = Int(du_dn*conj(v), t)*l  # int -du/dn *conj(v) dl
        I = -1j*k*d2*I_uv - I_duv
        return I

    def RHS(self) -> complex:
        raise NotImplementedError
        
class NtDNon_Local_waveguide:
    def __init__(self, R: float, d_2: float, n: int, N_MODES: int):
        self.R = R
        self.mode_n = n
        self.d_2 = d_2
        self.N_MODES = N_MODES
    
    def LHS(self, edge_u, edge_v, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        d2 = self.d_2
        d_n = d_phi
        d_m = d_psi
        R = self.R
        N_MODES = self.N_MODES

        P = edge_u["P"]

        t = np.linspace(0, 1, N_POINTS)
        N = edge_u["N"]
        T = edge_u["T"]
        l = edge_u["l"]
        x = P + np.outer(t,T)*l
        
        P_v = edge["P"]

        t = np.linspace(0, 1, N_POINTS)
        N = edge["N"]
        T = edge["T"]
        l = edge["l"]
        x = P + np.outer(t,T)*l

        theta = np.linspace(0, 2*np.pi, N_POINTS, endpoint=False)
        u_r = np.column_stack([np.cos(theta), np.sin(theta)])
        N = u_r
        x = R*u_r

        # u = np.where(theta_1_u < theta < theta_2_u, exp(1j*k*dot(x, d_n)), 0)
        mask_u = (theta_1_u <= theta) & (theta <= theta_2_u)
        u = np.zeros_like(theta, dtype=np.complex128)
        u[mask_u] = exp(1j*k*dot(x[mask_u, :], d_n))
        # v = np.where(theta_1_v < theta < theta_2_v, exp(1j*k*dot(x, d_m)), 0)
        mask_v = (theta_1_v <= theta) & (theta <= theta_2_v)
        v = np.zeros_like(theta, dtype=np.complex128)
        v[mask_v] = exp(1j*k*dot(x[mask_v, :], d_m))

        du_dn = 1j*k*dot(N, d_n)*u
        dv_dn = 1j*k*dot(N, d_m)*v

        NtD_du = NewmanntoDirichlet(y, du_dn, k, R, N_MODES)
        NtD_dv = NewmanntoDirichlet(y, dv_dn, k, R, N_MODES)
        
        I_Nudv = Int(NtD_du*conj(dv_dn), theta)*R
        I_NuNv = Int(NtD_du*conj(NtD_dv), theta)*R
        I_uNv  = Int(u*conj(NtD_dv), theta)*R
        I_Nuv  = Int(NtD_du*conj(v), theta)*R
        
        I = I_Nudv - d2*1j*k*(I_NuNv - I_Nuv - I_uNv)
        return I



def NewmanntoDirichlet(y, df_dy, k, H, M):

    dfn = np.zeros(M, dtype=np.complex128)
    dfn[0] = Int( df_dy*1/np.sqrt(H), y )
    for n in range(1,M):
        dfn[n] = Int( df_dy*cos(n*pi*y/H)/np.sqrt(H/2), y )
    
    f_y = 1/(1j*k)*dfn[0]/np.sqrt(H)*np.ones_like(y) + sum([ 1/(1j*np.sqrt(complex(k**2 - (n*pi/H)**2)))*dfn[n]*cos(n*pi*y/H)/np.sqrt(H/2) for n in range(1,M)])
    return f_y




# def num_RHS( k, P, Q, N, H, s, d_m, d2=0, Nt = 100, Np=15):
#     l = norm(Q-P)
#     t = np.linspace(0,1,Nt)
#     x = P + np.outer(t,Q-P)
#     psi_m = exp(1j*k*dot(x,d_m))
#     beta = sqrt(complex(k**2 - (s*pi/H)**2 ))
#     u_inc = exp(1j*beta*x[:,0])*cos(s*pi*x[:,1]/H)

#     grad_psi_m_N = 1j*k*dot(N,d_m)*exp(1j*k*dot(x,d_m))

#     N_gradpsi_m = NewmanntoDirichlet(x[:,1], grad_psi_m_N, k, H, Np)

#     I = -2*Int( u_inc*conj(grad_psi_m_N) - d2*1j*k*u_inc*conj(N_gradpsi_m-psi_m), t)*l
    
#     return I


# def test_RHS():
#     H=1
#     R= 10
#     P = np.array([-R,-H])
#     Q = np.array([-R,H])

#     l = norm(Q-P)
#     T = (Q - P)/l
#     N = np.array([0,-1])
#     M = (P+Q)/2
    

#     Edge = namedtuple('Edge',['P','Q','N','T', 'M', 'l'])
#     E = Edge(P,Q,N,T, M, l)

#     k = 8.
#     d_m = [1,1]
#     d_m = np.array(d_m)/norm(d_m)

#     TestFunction = namedtuple('TestFunction',['d','k'])
#     psi_m = TestFunction(d=d_m,k=k)

#     d2 = 0.5

#     t= 1

#     I_exact = exact_RHS_broken(psi_m, E, k, H, d2, t)
#     I_num = num_RHS( k, P, Q, N, H, t, d_m, d2=d2, Nt=N_points)
#     assert np.isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'


