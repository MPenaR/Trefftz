from cases.open_space.kernels import NtDLocal_circle, I_uv, I_duv, Fdudn, Fu, I_Nuv, I_uNv, I_NuNv
from cases.open_space.numerical_kernels import numerical_I_uv, numerical_I_duv, num_Fdudn, num_Fu, num_I_Nuv, num_I_uNv, num_I_NuNv
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


@pytest.mark.slow
@pytest.mark.parametrize('d_m', directions )
@pytest.mark.parametrize('d_n', directions )
def test_IuNv(d_m, d_n):
    
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

    I_exact = I_uNv(edge_u=E[0], edge_v=E[1], d_u=d_n, d_v=d_m, k=k, R=R, N_modes=60, NtD_modes=NtD_modes)
    I_num = num_I_uNv(edge_u=E[0], edge_v=E[1], d_u=d_n, d_v=d_m, k=k, R=R, NtD_modes=NtD_modes)
    assert isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'


@pytest.mark.slow
@pytest.mark.parametrize('d_m', directions )
@pytest.mark.parametrize('d_n', directions )
def test_IuNv(d_m, d_n):
    
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

    I_exact = I_NuNv(edge_u=E[0], edge_v=E[1], d_u=d_n, d_v=d_m, k=k, R=R, N_modes=60, NtD_modes=NtD_modes)
    I_num = num_I_NuNv(edge_u=E[0], edge_v=E[1], d_u=d_n, d_v=d_m, k=k, R=R, NtD_modes=NtD_modes)
    assert isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'


