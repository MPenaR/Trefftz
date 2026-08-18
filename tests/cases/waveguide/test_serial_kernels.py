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

# def num_NtDLocal_LHS(k, P, Q, N, H, d_n, d_m, d2=0, Nt = 100, Np=15) -> complex:
#     l = norm(Q-P)
#     t = np.linspace(0,1,Nt)
#     x = P + np.outer(t,Q-P)
#     phi_n = exp(1j*k*dot(x,d_n))
#     psi_m = exp(1j*k*dot(x,d_m))
#     grad_phi_n_N = 1j*k*dot(N,d_n)*exp(1j*k*dot(x,d_n))
#     # grad_psi_m_N = 1j*k*dot(N,d_m)*exp(1j*k*dot(x,d_m))

#     I = -Int(grad_phi_n_N*conj(psi_m) + d2*1j*k*phi_n*conj(psi_m), t)*l

#     return I


# @pytest.mark.parametrize('d_m', directions )
# @pytest.mark.parametrize('d_n', directions )
# def test_NtD_local_LHS(d_m, d_n):
#     P = array([0,1])
#     Q = array([3,1])
#     l = norm(P-Q)
#     T = (Q - P)/l
#     N = array([0,1])
#     M = (P + Q)/2
    
#     E = np.zeros((), dtype=edge_dtype)
#     E["P"] = P
#     E["Q"] = Q
#     E["N"] = N
#     E["T"] = T
#     E["M"] = M
#     E["l"] = l

#     k = 8.
#     d_n = array(d_n)/norm(d_n)
#     d_m = array(d_m)/norm(d_m)

#     H = 1.

#     d2 = 0.5
#     I_exact = NtDLocal(R=2, d_2=d2, n=1, H=H).LHS(edge=E, d_phi=d_n, d_psi=d_m, k=k)
#     I_num = numerical_NtDLocal(R=2, d_2=d2, n=1).LHS(edge=E, d_phi=d_n, d_psi=d_m, k=k)
    
#     # I_num = num_NtDLocal_LHS(k, P, Q, N, H, d_n, d_m, d2=d2,  Nt=N_POINTS)
#     assert isclose(I_num, I_exact, TOL, TOL), f'{I_exact=}, {I_num=}'
