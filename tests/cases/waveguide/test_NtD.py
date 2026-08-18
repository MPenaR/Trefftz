from ...data import DIRECTIONS, TOL, k
from .data import NTD_MODES, H, SEGMENT1_SIGMA_L, SEGMENT1_SIGMA_R
import cases.waveguide.NeumannToDirichlet as exact
import cases.waveguide.numerical_NtD as numerical
from numpy.testing import assert_allclose
import numpy as np
import pytest

N_POINTS = 10**5

@pytest.mark.parametrize('edge', (SEGMENT1_SIGMA_L, SEGMENT1_SIGMA_R))
@pytest.mark.parametrize('d', DIRECTIONS)
@pytest.mark.parametrize('t', range(NTD_MODES))
def test_ntd(edge, d, t):
    I_exact = exact.ntd(edge=edge, d=d, k=k, H=H, t=t)
    x = edge["P"][0]
    y = np.linspace(0, H, N_POINTS)
    r = np.column_stack(( x*np.ones_like(y), y))
    u = np.exp(1j*k*np.dot(r, d))
    N = edge["N"]
    mask = ( y >= edge["P"][1] ) & ( y <= edge["Q"][1] )
    du_dn = np.where(mask, 1j*k*np.dot(d, N)*u, 0.)

    I_num = numerical.ntd(y, du_dn, k=k, H=H, t=t)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)


@pytest.mark.parametrize('edge', (SEGMENT1_SIGMA_L, SEGMENT1_SIGMA_R))
@pytest.mark.parametrize('t', range(NTD_MODES))
def test_block_ntd(edge, t):
    I_exact = exact.ntd(edge=edge, d=DIRECTIONS, k=k, H=H, t=t)
    
    x = edge["P"][0]
    y = np.linspace(0, H, N_POINTS)
    r = np.column_stack(( x*np.ones_like(y), y))
    N = edge["N"]
    mask = ( y >= edge["P"][1] ) & ( y <= edge["Q"][1] )
    
    I_num = np.zeros_like(I_exact)
    for n, d in enumerate(DIRECTIONS):
        u = np.exp(1j*k*np.dot(r, d))
        du_dn = np.where(mask, 1j*k*np.dot(d, N)*u, 0.)
        I_num[n] = numerical.ntd(y, du_dn, k=k, H=H, t=t)
    
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)
