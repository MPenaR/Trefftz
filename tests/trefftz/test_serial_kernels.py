from trefftz.dg.serial_kernels import NeumannKernel, UltraWeakKernel, SIGN, I_uv_arc, I_uv, I_duv_arc
from trefftz.dg.numerical_kernels import numNeumannKernel, numUltraWeakKernel, numI_uv_arc, numI_uv, numI_duv_arc
import pytest
from numpy import linspace, sin, cos, pi, isclose, array
from numpy.linalg import norm
import numpy as np
from trefftz.mesh.core2 import edge_dtype


TOL = 1E-5
N_POINTS = int(1E5)


NTH = 3
directions = [(cos(th), sin(th)) for th in linspace(0, pi/2, NTH, endpoint=False)]


@pytest.mark.parametrize('d_m', directions )
@pytest.mark.parametrize('d_n', directions )
def test_Iuv(d_m, d_n):
    P = array([3,3])
    Q = array([1,1])
    l = norm(Q-P)
    T = (Q - P)/l
    N = array([-T[1], T[0]])
    M = (P + Q)/2    
    E = np.zeros((), dtype=edge_dtype)
    E["P"] = P
    E["Q"] = Q
    E["N"] = N
    E["T"] = T
    E["M"] = M
    E["l"] = l


    k = 8.
    d_n = array(d_n)
    d_m = array(d_m)

    I_exact = I_uv(edge=E, d_u=d_n, d_v=d_m, k=k)
    I_num = numI_uv(edge=E, d_u=d_n, d_v=d_m, k=k)
    assert isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'


@pytest.mark.parametrize('d_m', directions )
@pytest.mark.parametrize('d_n', directions )
def test_Iuv_arc(d_m, d_n):
    
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

    I_exact = I_uv_arc(edge=E, d_u=d_n, d_v=d_m, k=k, R=R)
    I_num = numI_uv_arc(edge=E, d_u=d_n, d_v=d_m, k=k, R=R)
    assert isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'


@pytest.mark.parametrize('d_m', directions )
@pytest.mark.parametrize('d_n', directions )
def test_Iduv_arc(d_m, d_n):
    
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

    I_exact = I_duv_arc(edge=E, d_u=d_n, d_v=d_m, k=k, R=R)
    I_num = numI_duv_arc(edge=E, d_u=d_n, d_v=d_m, k=k, R=R)
    assert isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'




@pytest.mark.parametrize('d_m', directions )
@pytest.mark.parametrize('d_n', directions )
def test_UltraWeakKernel(d_m,d_n):
    P = array([3,3])
    Q = array([1,1])
    l = norm(Q-P)
    T = (Q - P)/l
    N = array([-T[1], T[0]])
    M = (P + Q)/2    
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

    a = 0.5
    b = 0.5

    I_exact = UltraWeakKernel(a=a, b=b).LHS(edge=E, d_phi=d_n, d_psi=d_m, k=k, sign=SIGN.PP)
    I_num = numUltraWeakKernel(a=a, b=b).LHS(edge=E, d_phi=d_n, d_psi=d_m, k=k)
    assert isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'


@pytest.mark.parametrize('d_m', directions )
@pytest.mark.parametrize('d_n', directions )
def test_NeumannLHS(d_m, d_n):
    P = array([0,1])
    Q = array([3,1])
    l = norm(P-Q)
    T = (Q - P)/l
    N = array([0,1])
    M = (P + Q)/2
    
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

    d1 = 0.5
    I_exact = NeumannKernel(d_1=d1).LHS(edge=E, d_phi=d_n, d_psi=d_m, k=k)
    I_num = numNeumannKernel(d_1=d1).LHS(edge=E, d_phi=d_n, d_psi=d_m, k=k)
    
    assert isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'