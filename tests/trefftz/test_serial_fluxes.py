import pytest
from numpy import linspace, isclose, array
from numpy.linalg import norm
import numpy as np
from trefftz.mesh.core2 import edge_dtype
import trefftz.dg.serial_fluxes as exact
import trefftz.dg.numerical_fluxes as numerical



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