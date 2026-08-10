import pytest
from numpy.testing import assert_allclose
import trefftz.dg.serial_fluxes as exact
import trefftz.dg.numerical_fluxes as numerical
from ..data import k, DIRECTIONS, EDGE_1, ARC_1, TOL


@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_UltraWeak(d_m, d_n):
    a = 0.5
    b = 0.5
    I_exact = exact.UltraWeakFlux(a=a, b=b).LHS(edge=EDGE_1, d_phi=d_n, d_psi=d_m, k=k, sign=exact.SIGN.PP)
    I_num = numerical.UltraWeakFlux(a=a, b=b).LHS(edge=EDGE_1, d_phi=d_n, d_psi=d_m, k=k, sign=exact.SIGN.PP)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)  # sign is not really checked as both numerical and exact are using the "exact" understanding of sign


@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_NeumannLHS(d_m, d_n):
    d1 = 0.5
    I_exact = exact.NeumannFlux(d_1=d1).LHS(edge=EDGE_1, d_phi=d_n, d_psi=d_m, k=k)
    I_num = numerical.NeumannFlux(d_1=d1).LHS(edge=EDGE_1, d_phi=d_n, d_psi=d_m, k=k)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)


@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_DirichletLHS(d_m, d_n):
    a = 0.5
    I_exact = exact.DirichletFlux(a=a).LHS(edge=EDGE_1, d_phi=d_n, d_psi=d_m, k=k)
    I_num = numerical.DirichletFlux(a=a).LHS(edge=EDGE_1, d_phi=d_n, d_psi=d_m, k=k)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)


@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_inc', DIRECTIONS)
def test_DirichletRHS(d_m, d_inc):
    a = 0.5
    data = {"d_inc": d_inc}
    I_exact = exact.DirichletFlux(a=a, data=data).RHS(edge=EDGE_1, d_psi=d_m, k=k)
    I_num = numerical.DirichletFlux(a=a, data=data).RHS(edge=EDGE_1, d_psi=d_m, k=k)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)


@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_CircularDirichletLHS(d_m, d_n):
    a = 0.5
    I_exact = exact.DirichletFlux(a=a).LHS(edge=ARC_1, d_phi=d_n, d_psi=d_m, k=k)
    I_num = numerical.CircularDirichletFlux(a=a).LHS(edge=ARC_1, d_phi=d_n, d_psi=d_m, k=k)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)


@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_inc', DIRECTIONS)
def test_CircularDirichletRHS(d_m, d_inc):
    a = 0.5
    data = {"d_inc": d_inc}
    I_exact = exact.DirichletFlux(a=a, data=data).RHS(edge=ARC_1, d_psi=d_m, k=k)
    I_num = numerical.CircularDirichletFlux(a=a, data=data).RHS(edge=ARC_1, d_psi=d_m, k=k)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)
