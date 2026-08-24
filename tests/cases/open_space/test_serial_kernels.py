import cases.open_space.serial_kernels as exact
import cases.open_space.numerical.kernels as numerical
import pytest
from numpy.testing import assert_allclose
from ...data import TOL, DIRECTIONS, k
from .data import ARC1_SIGMA, ARC2_SIGMA, NTD_MODES, JACOBI_ANGER_MODES

@pytest.mark.slow
@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_IuNv(d_m, d_n):
    arc_u = ARC1_SIGMA
    arc_v = ARC2_SIGMA
    I_exact = exact.I_uNv(arc_u=arc_u, arc_v=arc_v, d_u=d_n, d_v=d_m, k=k, NtD_modes= NTD_MODES, JA_modes = JACOBI_ANGER_MODES)
    I_num = numerical.I_uNv(arc_u=arc_u, arc_v=arc_v, d_u=d_n, d_v=d_m, k=k, NtD_modes= NTD_MODES)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)

@pytest.mark.slow
@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_INuv(d_m, d_n):
    arc_u = ARC1_SIGMA
    arc_v = ARC2_SIGMA
    I_exact = exact.I_Nuv(arc_u=arc_u, arc_v=arc_v, d_u=d_n, d_v=d_m, k=k, NtD_modes= NTD_MODES, JA_modes = JACOBI_ANGER_MODES)
    I_num = numerical.I_Nuv(arc_u=arc_u, arc_v=arc_v, d_u=d_n, d_v=d_m, k=k, NtD_modes= NTD_MODES)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)


@pytest.mark.slow
@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_INuNv(d_m, d_n):
    arc_u = ARC1_SIGMA
    arc_v = ARC2_SIGMA
    I_exact = exact.I_NuNv(arc_u=arc_u, arc_v=arc_v, d_u=d_n, d_v=d_m, k=k, NtD_modes= NTD_MODES, JA_modes = JACOBI_ANGER_MODES)
    I_num = numerical.I_NuNv(arc_u=arc_u, arc_v=arc_v, d_u=d_n, d_v=d_m, k=k, NtD_modes= NTD_MODES)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)


@pytest.mark.slow
@pytest.mark.parametrize('d_m', DIRECTIONS)
@pytest.mark.parametrize('d_n', DIRECTIONS)
def test_INuNv(d_m, d_n):
    arc_u = ARC1_SIGMA
    arc_v = ARC2_SIGMA
    I_exact = exact.I_NuNv(arc_u=arc_u, arc_v=arc_v, d_u=d_n, d_v=d_m, k=k, NtD_modes= NTD_MODES, JA_modes = JACOBI_ANGER_MODES)
    I_num = numerical.I_NuNv(arc_u=arc_u, arc_v=arc_v, d_u=d_n, d_v=d_m, k=k, NtD_modes= NTD_MODES)
    assert_allclose(I_num, I_exact, rtol=TOL, atol=TOL)
