from ...data import DIRECTIONS, TOL, k
from .data import NTD_MODES,  ARC1_SIGMA, JACOBI_ANGER_MODES
import cases.open_space.NeumannToDirichlet as exact
import cases.open_space.numerical_NtD as numerical
from numpy.testing import assert_allclose
import numpy as np
from numpy import exp, dot
import pytest

N_POINTS = 10**7

@pytest.mark.parametrize('d', DIRECTIONS)
@pytest.mark.parametrize('t', range(NTD_MODES))
def test_ntd(d, t):
    I_exact = exact.ntd(arc=ARC1_SIGMA, d=d, k=k, t=t, JA_modes=JACOBI_ANGER_MODES)
    theta_1 = ARC1_SIGMA["theta_1"]
    theta_2 = ARC1_SIGMA["theta_2"]

    R = ARC1_SIGMA["R"]

    theta = np.linspace(0, 2*np.pi, N_POINTS, endpoint=False)
    u_r = np.column_stack([np.cos(theta), np.sin(theta)])
    N = u_r
    x = R*u_r

    mask = (theta_1 <= theta) & (theta <= theta_2)
    u = np.where(mask, exp(1j*k*dot(x, d)), 0.)
    du_dn = 1j*k*dot(N, d)*u
    I_num = numerical.ntd(theta, du_dn, k=k, R=R, t=t)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)


@pytest.mark.parametrize('t', range(NTD_MODES))
def test_ntd_block(t):
    I_exact = exact.ntd(arc=ARC1_SIGMA, d=DIRECTIONS, k=k, t=t, JA_modes=JACOBI_ANGER_MODES)
    I_num = np.zeros_like(I_exact)
    theta_1 = ARC1_SIGMA["theta_1"]
    theta_2 = ARC1_SIGMA["theta_2"]

    R = ARC1_SIGMA["R"]

    theta = np.linspace(0, 2*np.pi, N_POINTS, endpoint=False)
    u_r = np.column_stack([np.cos(theta), np.sin(theta)])
    N = u_r
    x = R*u_r

    mask = (theta_1 <= theta) & (theta <= theta_2)
    for n, d in enumerate(DIRECTIONS):
        u = np.where(mask, exp(1j*k*dot(x, d)), 0.)
        du_dn = 1j*k*dot(N, d)*u
        I_num[n] = numerical.ntd(theta, du_dn, k=k, R=R, t=t)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)