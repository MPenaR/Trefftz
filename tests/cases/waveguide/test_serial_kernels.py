import cases.waveguide.serial_kernels as exact
import cases.waveguide.numerical_kernels as numerical
import pytest
from numpy.testing import assert_allclose
from ...data import TOL, DIRECTIONS, k
from .data import EDGE1_SIGMA_L,EDGE2_SIGMA_L, EDGE1_SIGMA_R, EDGE2_SIGMA_R, H, NTD_MODES

@pytest.mark.parametrize('edges', [(EDGE1_SIGMA_L, EDGE2_SIGMA_L),
                                   (EDGE1_SIGMA_R, EDGE2_SIGMA_R)])
@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_IuNv(edges, d_m, d_n):
    edge_u, edge_v = edges
    I_exact = exact.I_uNv(segment_u=edge_u, segment_v=edge_v, d_u=d_n, d_v=d_m, k=k, H=H, NtD_modes= NTD_MODES)
    I_num = numerical.I_uNv(segment_u=edge_u, segment_v=edge_v, d_u=d_n, d_v=d_m, k=k, H=H, NtD_modes= NTD_MODES)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)


@pytest.mark.parametrize('edges', [(EDGE1_SIGMA_L, EDGE2_SIGMA_L),
                                   (EDGE1_SIGMA_R, EDGE2_SIGMA_R)])
@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_INuv(edges, d_m, d_n):
    edge_u, edge_v = edges
    I_exact = exact.I_Nuv(segment_u=edge_u, segment_v=edge_v, d_u=d_n, d_v=d_m, k=k, H=H, NtD_modes= NTD_MODES)
    I_num = numerical.I_Nuv(segment_u=edge_u, segment_v=edge_v, d_u=d_n, d_v=d_m, k=k, H=H, NtD_modes= NTD_MODES)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)


@pytest.mark.parametrize('edges', [(EDGE1_SIGMA_L, EDGE2_SIGMA_L),
                                   (EDGE1_SIGMA_R, EDGE2_SIGMA_R)])
@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_INuNv(edges, d_m, d_n):
    edge_u, edge_v = edges
    I_exact = exact.I_NuNv(segment_u=edge_u, segment_v=edge_v, d_u=d_n, d_v=d_m, k=k, H=H, NtD_modes= NTD_MODES)
    I_num = numerical.I_NuNv(segment_u=edge_u, segment_v=edge_v, d_u=d_n, d_v=d_m, k=k, H=H, NtD_modes= NTD_MODES)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)


@pytest.mark.parametrize('edges', [(EDGE1_SIGMA_L, EDGE2_SIGMA_L),
                                   (EDGE1_SIGMA_R, EDGE2_SIGMA_R)])
@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_INudv(edges, d_m, d_n):
    edge_u, edge_v = edges
    I_exact = exact.I_Nudv(segment_u=edge_u, segment_v=edge_v, d_u=d_n, d_v=d_m, k=k, H=H, NtD_modes= NTD_MODES)
    I_num = numerical.I_Nudv(segment_u=edge_u, segment_v=edge_v, d_u=d_n, d_v=d_m, k=k, H=H, NtD_modes= NTD_MODES)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)