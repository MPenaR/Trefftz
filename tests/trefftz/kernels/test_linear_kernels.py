import pytest
from numpy.testing import assert_allclose
import trefftz.dg.kernels.serial.linear_kernels as exact
import trefftz.dg.numerical.kernels as numerical
from ...data import TOL, DIRECTIONS, EDGE_1, k, n_1, n_2, planewave


@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_Iuv_segment(d_m, d_n):
    I_exact = exact.I_uv(segment=EDGE_1, d_u=d_n, d_v=d_m, k=k, n_u=n_1, n_v=n_2)
    I_num = numerical.I_uv_segment(segment=EDGE_1, d_u=d_n, d_v=d_m, k=k, n_u=n_1, n_v=n_2)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)

@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_Iduv_segment(d_m, d_n):
    I_exact = exact.I_duv(segment=EDGE_1, d_u=d_n, d_v=d_m, k=k, n_u=n_1, n_v=n_2)
    I_num = numerical.I_duv_segment(segment=EDGE_1, d_u=d_n, d_v=d_m, k=k, n_u=n_1, n_v=n_2)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)

@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_Iudv_segment(d_m, d_n):
    I_exact = exact.I_udv(segment=EDGE_1, d_u=d_n, d_v=d_m, k=k, n_u=n_1, n_v=n_2)
    I_num = numerical.I_udv_segment(segment=EDGE_1, d_u=d_n, d_v=d_m, k=k, n_u=n_1, n_v=n_2)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)

@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_Idudv_segment(d_m, d_n):
    I_exact = exact.I_dudv(segment=EDGE_1, d_u=d_n, d_v=d_m, k=k, n_u=n_1, n_v=n_2)
    I_num = numerical.I_dudv_segment(segment=EDGE_1, d_u=d_n, d_v=d_m, k=k, n_u=n_1, n_v=n_2)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)

@pytest.mark.parametrize('d_m', DIRECTIONS)
def test_I_pw_v_segment(d_m):
    I_exact = exact.I_pw_v(segment_v=EDGE_1, plane_wave=planewave, d_v=d_m, k=k, n_v=n_2)
    I_num = numerical.I_pw_v_segment(segment=EDGE_1, planewave=planewave, d_v=d_m, k=k, n_v=n_2)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)

@pytest.mark.parametrize('d_m', DIRECTIONS)
def test_I_pw_dv_segment(d_m):
    I_exact = exact.I_pw_dv(segment_v=EDGE_1, plane_wave=planewave, d_v=d_m, k=k, n_v=n_2)
    I_num = numerical.I_pw_dv_segment(segment=EDGE_1, planewave=planewave, d_v=d_m, k=k, n_v=n_2)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)

@pytest.mark.parametrize('d_m', DIRECTIONS)
def test_I_dpw_v_segment(d_m):
    I_exact = exact.I_dpw_v(segment_v=EDGE_1, plane_wave=planewave, d_v=d_m, k=k, n_v=n_2)
    I_num = numerical.I_dpw_v_segment(segment=EDGE_1, planewave=planewave, d_v=d_m, k=k, n_v=n_2)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)