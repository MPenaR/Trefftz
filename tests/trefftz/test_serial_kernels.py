import pytest
from numpy import isclose, array
from numpy.linalg import norm
import numpy as np
from trefftz.mesh.core2 import edge_dtype
import trefftz.dg.serial_kernels as exact
import trefftz.dg.numerical_kernels as numerical
from tests.parameters import TOL, DIRECTIONS


@pytest.mark.parametrize('d_m', DIRECTIONS )
@pytest.mark.parametrize('d_n', DIRECTIONS )
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

    I_exact = exact.I_uv(edge=E, d_u=d_n, d_v=d_m, k=k)
    I_num = numerical.I_uv(edge=E, d_u=d_n, d_v=d_m, k=k)
    assert isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'


@pytest.mark.parametrize('d_m', DIRECTIONS )
@pytest.mark.parametrize('d_n', DIRECTIONS )
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

    I_exact = exact.I_uv_arc(edge=E, d_u=d_n, d_v=d_m, k=k)
    I_num = numerical.I_uv_arc(edge=E, d_u=d_n, d_v=d_m, k=k, R=R)
    assert isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'


@pytest.mark.parametrize('d_m', DIRECTIONS )
@pytest.mark.parametrize('d_n', DIRECTIONS )
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

    I_exact = exact.I_duv_arc(edge=E, d_u=d_n, d_v=d_m, k=k)
    I_num = numerical.I_duv_arc(edge=E, d_u=d_n, d_v=d_m, k=k, R=R)
    assert isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'




