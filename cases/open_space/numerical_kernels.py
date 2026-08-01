from numpy import linspace, outer, sin, cos, pi, exp, dot, conj, isclose, array
from numpy.lib.scimath import sqrt
from numpy.linalg import norm
from numpy import trapezoid as Int
import numpy as np
from scipy.special import hankel1, h1vp
from trefftz.numpy_types import float_array, complex_array

N_POINTS = int(1E5)
NtD_MODES = 20



class NtDLocal_circle:
    def __init__(self, R: float, d_2: float, n: int):
        self.R = R
        self.mode_n = n
        self.d_2 = d_2

    
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        d2 = self.d_2
        d_n = d_phi
        d_m = d_psi
        R = self.R

        P = edge["P"]
        Q = edge["Q"]
        theta_1 = np.atan2(P[1], P[0])
        theta_2 = np.atan2(Q[1], Q[0])

        theta = np.linspace(theta_1, theta_2, N_POINTS)
        u_r = np.column_stack([np.cos(theta), np.sin(theta)])
        N = u_r
        x = R*u_r
        u = exp(1j*k*dot(x, d_n))
        v = exp(1j*k*dot(x, d_m))
        du_dn = 1j*k*dot(N, d_n)*u
        I_uv = -1j*k*d2*Int(u*conj(v), theta)*R  # int -i*k*d_2*u*conj(v) dl 
        I_duv = -Int(du_dn*conj(v), theta)*R  # int -du/dn *conj(v) dl
        I = I_uv + I_duv
        return I

    def RHS(self) -> complex:
        raise NotImplementedError
        

def numerical_I_uv(edge, d_u: float_array, d_v: float_array, k: float, R: float) -> complex:
    P = edge["P"]
    Q = edge["Q"]
    theta_1 = np.atan2(P[1], P[0])
    theta_2 = np.atan2(Q[1], Q[0])

    theta = np.linspace(theta_1, theta_2, N_POINTS)
    u_r = np.column_stack([np.cos(theta), np.sin(theta)])
    x = R*u_r
    u = exp(1j*k*dot(x, d_u))
    v = exp(1j*k*dot(x, d_v))
    I = Int(u*conj(v), theta)*R  # int u*conj(v) dl 
    return I

def numerical_I_duv(edge, d_u: float_array, d_v: float_array, k: float, R: float) -> complex:
    P = edge["P"]
    Q = edge["Q"]
    theta_1 = np.atan2(P[1], P[0])
    theta_2 = np.atan2(Q[1], Q[0])

    theta = np.linspace(theta_1, theta_2, N_POINTS)
    u_r = np.column_stack([np.cos(theta), np.sin(theta)])
    N = u_r
    x = R*u_r

    u = exp(1j*k*dot(x, d_u))
    v = exp(1j*k*dot(x, d_v))

    du_dn = 1j*k*dot(N, d_u)*u

    I = Int(du_dn*conj(v), theta)*R  # int du/dn*conj(v) dl 
    return I


class NtDNon_Local_circle:
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

        P_u = edge_u["P"]
        Q_u = edge_u["Q"]
        theta_1_u = np.atan2(P_u[1], P_u[0])
        theta_2_u = np.atan2(Q_u[1], Q_u[0])
        
        P_v = edge_v["P"]
        Q_v = edge_v["Q"]
        theta_1_v = np.atan2(P_v[1], P_v[0])
        theta_2_v = np.atan2(Q_v[1], Q_v[0])

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

        NtD_du = NewmanntoDirichlet(theta, du_dn, k, R, N_MODES)
        NtD_dv = NewmanntoDirichlet(theta, dv_dn, k, R, N_MODES)
        
        I_Nudv = Int(NtD_du*conj(dv_dn), theta)*R
        I_NuNv = Int(NtD_du*conj(NtD_dv), theta)*R
        I_uNv  = Int(u*conj(NtD_dv), theta)*R
        I_Nuv  = Int(NtD_du*conj(v), theta)*R
        
        I = I_Nudv - d2*1j*k*(I_NuNv - I_Nuv - I_uNv)
        return I




def FourierCoefficient(theta: float_array, f: float_array, t: int) -> complex:
    r'''Computes the t Fourier coefficient defined as:
    
    .. math ::
        \frac{1}{2\pi}*\int_0^{2\pi}f(\theta)e^{-it\theta}\,\mathrm{d}\theta
    '''
    return 1/sqrt(2*np.pi)*Int(f*exp(-1j*t*theta), theta)


def num_Fdudn(edge, d_u: float_array, k: float, R: float, t: int) -> complex:
    r'''Computes the t Fourier coefficient of 
    .. math ::
    \nabla u\cdot \mathbf{n} 
    supported at a single edge.
    '''
    P = edge["P"]
    Q = edge["Q"]
    theta_1 = np.atan2(P[1], P[0])
    theta_2 = np.atan2(Q[1], Q[0])

    theta = np.linspace(theta_1, theta_2, N_POINTS)
    u_r = np.column_stack([np.cos(theta), np.sin(theta)])
    N = u_r
    x = R*u_r

    u = exp(1j*k*dot(x, d_u))
    du_dn = 1j*k*dot(N, d_u)*u

    I = FourierCoefficient(theta, du_dn, t)

    return I

def num_Fu(edge, d: float_array, k: float, R: float, t: int) -> complex:
    r'''Computes the t Fourier coefficient of 
    .. math ::
    \nabla u\cdot \mathbf{n} 
    supported at a single edge.
    '''
    P = edge["P"]
    Q = edge["Q"]
    theta_1 = np.atan2(P[1], P[0])
    theta_2 = np.atan2(Q[1], Q[0])

    theta = np.linspace(theta_1, theta_2, N_POINTS)
    u_r = np.column_stack([np.cos(theta), np.sin(theta)])
    x = R*u_r

    u = exp(1j*k*dot(x, d))
    I = FourierCoefficient(theta, u, t)

    return I

    # return 1/(2*np.pi)*R*Int(du_dn*exp(-1j*t*theta), theta)

def NewmanntoDirichlet(theta: float_array, df_dn: complex_array, k: float, R: float, M: int) -> complex_array:
    L = 2*np.pi*R
    # dfn = np.zeros(M, dtype=np.complex128)
    # dfn[0] = R*Int( df_dn*1/np.sqrt(L), theta )
    # for n in range(1,M):
    #     dfn[n] = Int( df_dy*cos(n*pi*y/H)/np.sqrt(H/2), y )

    # f_y = 1/(1j*k)*dfn[0]/np.sqrt(H)*np.ones_like(y) + sum([ 1/(1j*np.sqrt(complex(k**2 - (n*pi/H)**2)))*dfn[n]*cos(n*pi*y/H)/np.sqrt(H/2) for n in range(1,M)])

    dfn = np.zeros(2*M+1, dtype=np.complex128)
    for n in range(-M, M+1):
        i = n+M
        dfn[i] = R*Int(df_dn*exp(-1j*n*theta)/np.sqrt(L), theta)

    f = 1/k*sum((hankel1(n, k*R)/h1vp(n, k*R)*dfn[n+M]*exp(1j*n*theta)/np.sqrt(L) for n in range(-M, M+1)))

    return f






# def num_Radiating( k, P, Q, N, H, d_n, d_m, d2=0, Nt = 100, N_modes=15):
#     l = norm(Q-P)
#     t = np.linspace(0,1,Nt)
#     x = P + np.outer(t,Q-P)
#     phi_n = exp(1j*k*dot(x,d_n))
#     psi_m = exp(1j*k*dot(x,d_m))
#     grad_phi_n_N = 1j*k*dot(N,d_n)*exp(1j*k*dot(x,d_n))
#     grad_psi_m_N = 1j*k*dot(N,d_m)*exp(1j*k*dot(x,d_m))

#     N_gradphi_n = NewmanntoDirichlet(x[:,1], grad_phi_n_N, k, H, N_modes)
#     N_gradpsi_m = NewmanntoDirichlet(x[:,1], grad_psi_m_N, k, H, N_modes)

#     I = Int( N_gradphi_n*conj(grad_psi_m_N) - grad_phi_n_N*conj(psi_m), t)*l
#     I+= -d2*1j*k*Int((N_gradphi_n - phi_n)*conj(N_gradpsi_m - psi_m), t)*l
    
#     return I

# #@pytest.mark.xfail(reason="mixed up dimensions of the waveguide")
# @pytest.mark.parametrize(('d_m', 'd_n'), directions )
# def test_Radiating(d_m,d_n):
#     H=1
#     R= 10
#     P = np.array([R,0])
#     Q = np.array([R,H])

#     l = norm(Q-P)
#     T = (Q - P)/l
#     N = np.array([1,0])
#     M = (P+Q)/2

#     E = Edge(P,Q,N,T,M,l)

#     k = 8.
#     d_n = np.array(d_n)/norm(d_n)
#     d_m = np.array(d_m)/norm(d_m)

#     phi = Function(d=d_n,n=1)
#     psi = Function(d=d_m,n=1)
#     d_2 = 0.5

#     N_modes = 15
#     I_exact_local = Radiating_local(phi, psi, k, E, d_2)
#     I_exact_nonlocal = Radiating_nonlocal(phi=phi, psi=psi, k=k, edge_u=E, edge_v=E, d_2=d_2, N_modes=N_modes, H=H)
#     I_exact = I_exact_nonlocal + I_exact_local
#     I_num = num_Radiating( k, P, Q, N, H, d_n, d_m, d2=d_2,  Nt=N_POINTS, N_modes=N_modes)
#     assert np.isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'



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


# def test_RHS_broken():
#     H=1
#     R= 10
#     P = np.array([-R,-H])
#     Q = np.array([-R,H/3])

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
