import pytest
from numpy.testing import assert_allclose
import trefftz.dg.serial_kernels as exact
import trefftz.dg.numerical_kernels as numerical
from ..data import TOL, DIRECTIONS, EDGE_1, ARC_1, k


@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_Iuv(d_m, d_n):
    I_exact = exact.I_uv_segment(segment=EDGE_1, d_u=d_n, d_v=d_m, k=k)
    I_num = numerical.I_uv_segment(segment=EDGE_1, d_u=d_n, d_v=d_m, k=k)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)


@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_Iuv_arc(d_m, d_n):

    I_exact = exact.I_uv_arc(arc=ARC_1, d_u=d_n, d_v=d_m, k=k)
    I_num = numerical.I_uv_arc(arc=ARC_1, d_u=d_n, d_v=d_m, k=k)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)


@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_Iduv_arc(d_m, d_n):
    I_exact = exact.I_duv_arc(arc=ARC_1, d_u=d_n, d_v=d_m, k=k)
    I_num = numerical.I_duv_arc(arc=ARC_1, d_u=d_n, d_v=d_m, k=k)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)