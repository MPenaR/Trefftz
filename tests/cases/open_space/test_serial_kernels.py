from cases.open_space.kernels import NtDLocal_circle, I_uv, I_duv, Fdudn, Fu, I_Nuv
from cases.open_space.numerical_kernels import numerical_I_uv, numerical_I_duv, num_Fdudn, num_Fu, num_I_Nuv
import pytest
from numpy import linspace, outer, sin, cos, pi, exp, dot, conj, isclose, array
from numpy.lib.scimath import sqrt
from numpy.linalg import norm
from numpy import trapezoid as Int
import numpy as np
from trefftz.mesh.core2 import edge_dtype

TOL = 1E-7
N_POINTS = int(1E5)


NTH = 3
directions = [(cos(th), sin(th)) for th in linspace(0, pi/2, NTH, endpoint=False)]

@pytest.mark.parametrize('d_m', directions )
@pytest.mark.parametrize('d_n', directions )
def test_Iuv(d_m, d_n):
    
    theta_1 = np.pi*30/180
    theta_2 = np.pi*45/180
    R = 2.
    

    P = R*array([np.cos(theta_1), np.sin(theta_1)])
    Q = R*array([np.cos(theta_2), np.sin(theta_2)])
    l = norm(P-Q)
    T = (Q - P)/l
    N = array([0,1]) # meaningless, is a curved edge
    M = (P + Q)/2 # meaningless, is a curved edge
    
    E = np.zeros((), dtype=edge_dtype)
    E["P"] = P
    E["Q"] = Q
    E["N"] = N
    E["T"] = T
    E["M"] = M
    E["l"] = l

    k = 8.
    d_n = array(d_n)/norm(d_n)
    d_m = array(d_m)/norm(d_m)

    I_exact = I_uv(edge=E, d_u=d_n, d_v=d_m, k=k, R=R, N_modes=60)
    I_num = numerical_I_uv(edge=E, d_u=d_n, d_v=d_m, k=k, R=R)
    assert isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'


@pytest.mark.parametrize('d_m', directions )
@pytest.mark.parametrize('d_n', directions )
def test_Iduv(d_m, d_n):
    
    theta_1 = np.pi*30/180
    theta_2 = np.pi*45/180
    R = 2.
    

    P = R*array([np.cos(theta_1), np.sin(theta_1)])
    Q = R*array([np.cos(theta_2), np.sin(theta_2)])
    l = norm(P-Q)
    T = (Q - P)/l
    N = array([0,1]) # meaningless, is a curved edge
    M = (P + Q)/2 # meaningless, is a curved edge
    
    E = np.zeros((), dtype=edge_dtype)
    E["P"] = P
    E["Q"] = Q
    E["N"] = N
    E["T"] = T
    E["M"] = M
    E["l"] = l

    k = 8.
    d_n = array(d_n)/norm(d_n)
    d_m = array(d_m)/norm(d_m)

    I_exact = I_duv(edge=E, d_u=d_n, d_v=d_m, k=k, R=R, N_modes=40)
    I_num = numerical_I_duv(edge=E, d_u=d_n, d_v=d_m, k=k, R=R, )
    assert isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'


@pytest.mark.parametrize('t', [0, 1, 2] )
@pytest.mark.parametrize('d_n', directions )
def test_Fdudn(d_n, t):
    
    theta_1 = np.pi*30/180
    theta_2 = np.pi*45/180
    R = 2.
    

    P = R*array([np.cos(theta_1), np.sin(theta_1)])
    Q = R*array([np.cos(theta_2), np.sin(theta_2)])
    l = norm(P-Q)
    T = (Q - P)/l
    N = array([0,1]) # meaningless, is a curved edge
    M = (P + Q)/2 # meaningless, is a curved edge
    
    E = np.zeros((), dtype=edge_dtype)
    E["P"] = P
    E["Q"] = Q
    E["N"] = N
    E["T"] = T
    E["M"] = M
    E["l"] = l

    k = 8.
    d_n = array(d_n)/norm(d_n)

    I_exact = Fdudn(edge=E, d=d_n, k=k, R=R, N_modes=40, t=t)
    I_num = num_Fdudn(edge=E, d=d_n, k=k, R=R, t=t)
    assert isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'

@pytest.mark.parametrize('t', [0, 1, 2] )
@pytest.mark.parametrize('d_n', directions )
def test_Fu(d_n, t):
    
    theta_1 = np.pi*30/180
    theta_2 = np.pi*45/180
    R = 2.
    

    P = R*array([np.cos(theta_1), np.sin(theta_1)])
    Q = R*array([np.cos(theta_2), np.sin(theta_2)])
    l = norm(P-Q)
    T = (Q - P)/l
    N = array([0,1]) # meaningless, is a curved edge
    M = (P + Q)/2 # meaningless, is a curved edge
    
    E = np.zeros((), dtype=edge_dtype)
    E["P"] = P
    E["Q"] = Q
    E["N"] = N
    E["T"] = T
    E["M"] = M
    E["l"] = l

    k = 8.
    d_n = array(d_n)/norm(d_n)

    I_exact = Fu(edge=E, d=d_n, k=k, R=R, N_modes=40, t=t)
    I_num = num_Fu(edge=E, d=d_n, k=k, R=R, t=t)
    assert isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'


@pytest.mark.slow
@pytest.mark.parametrize('d_m', directions )
@pytest.mark.parametrize('d_n', directions )
def test_INuv(d_m, d_n):
    
    R = 3.

    thetas = [ (np.pi*30/180, np.pi*45/180 ),
               (np.pi*60/180, np.pi*90/180,)]

    E = np.zeros((2,), dtype=edge_dtype)

    for i, (theta_1, theta_2) in enumerate(thetas):
        P = R*array([np.cos(theta_1), np.sin(theta_1)])
        Q = R*array([np.cos(theta_2), np.sin(theta_2)])
        l = norm(P-Q)
        T = (Q - P)/l
        N = array([0,1]) # meaningless, is a curved edge
        M = (P + Q)/2 # meaningless, is a curved edge
        
        E[i]["P"] = P
        E[i]["Q"] = Q
        E[i]["N"] = N
        E[i]["T"] = T
        E[i]["M"] = M
        E[i]["l"] = l

    k = 8.
    d_n = array(d_n)/norm(d_n)
    d_m = array(d_m)/norm(d_m)

    NtD_modes = 3

    I_exact = I_Nuv(edge_u=E[0], edge_v=E[1], d_u=d_n, d_v=d_m, k=k, R=R, N_modes=60, NtD_modes=NtD_modes)
    I_num = num_I_Nuv(edge_u=E[0], edge_v=E[1], d_u=d_n, d_v=d_m, k=k, R=R, NtD_modes=NtD_modes)
    assert isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'




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
