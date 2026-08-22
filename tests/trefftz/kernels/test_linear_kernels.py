import pytest
from numpy.testing import assert_allclose
import trefftz.dg.kernels.serial.linear_kernels as exact
import trefftz.dg.numerical.kernels as numerical
from ...data import TOL, DIRECTIONS, EDGE_1,  k


@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_Iuv_segment(d_m, d_n):
    I_exact = exact.I_uv(segment=EDGE_1, d_u=d_n, d_v=d_m, k=k)
    I_num = numerical.I_uv_segment(segment=EDGE_1, d_u=d_n, d_v=d_m, k=k)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)

@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_Iduv_segment(d_m, d_n):
    I_exact = exact.I_duv(segment=EDGE_1, d_u=d_n, d_v=d_m, k=k)
    I_num = numerical.I_duv_segment(segment=EDGE_1, d_u=d_n, d_v=d_m, k=k)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)

@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_Iudv_segment(d_m, d_n):
    I_exact = exact.I_udv(segment=EDGE_1, d_u=d_n, d_v=d_m, k=k)
    I_num = numerical.I_udv_segment(segment=EDGE_1, d_u=d_n, d_v=d_m, k=k)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)