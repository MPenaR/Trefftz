import pytest
from numpy.testing import assert_allclose
import trefftz.dg.kernels.arc_kernels as exact
import trefftz.dg.numerical_kernels as numerical
from ..data import TOL, DIRECTIONS, ARC_1, k


@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_Iuv_arc(d_m, d_n):

    I_exact = exact.I_uv(arc=ARC_1, d_u=d_n, d_v=d_m, k=k)
    I_num = numerical.I_uv_arc(arc=ARC_1, d_u=d_n, d_v=d_m, k=k)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)


@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_Iduv_arc(d_m, d_n):
    I_exact = exact.I_duv(arc=ARC_1, d_u=d_n, d_v=d_m, k=k)
    I_num = numerical.I_duv_arc(arc=ARC_1, d_u=d_n, d_v=d_m, k=k)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)


@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_Iudv_arc(d_m, d_n):
    I_exact = exact.I_udv(arc=ARC_1, d_u=d_n, d_v=d_m, k=k)
    I_num = numerical.I_udv_arc(arc=ARC_1, d_u=d_n, d_v=d_m, k=k)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)