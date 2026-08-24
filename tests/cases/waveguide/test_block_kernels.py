import cases.waveguide.serial_kernels as serial
import cases.waveguide.block_kernels as block
import numpy as np
import pytest
from numpy.testing import assert_allclose
from ...data import TOL, DIRECTIONS, k
from .data import SEGMENT1_SIGMA_L, SEGMENT2_SIGMA_L, SEGMENT1_SIGMA_R, SEGMENT2_SIGMA_R, H, NTD_MODES


@pytest.mark.parametrize('edges', [(SEGMENT1_SIGMA_L, SEGMENT2_SIGMA_L),
                                   (SEGMENT1_SIGMA_R, SEGMENT2_SIGMA_R)])
def test_IuNv(edges):
    edge_u, edge_v = edges
    I_block = block.I_uNv(segment_u=edge_u, segment_v=edge_v, D_u=DIRECTIONS, D_v=DIRECTIONS, k=k, H=H, NtD_modes= NTD_MODES)
    I_serial = np.zeros_like(I_block)
    for n, d_n in enumerate(DIRECTIONS):
        for m, d_m in enumerate(DIRECTIONS):
            I_serial[m, n] = serial.I_uNv(segment_u=edge_u, segment_v=edge_v, d_u=d_n, d_v=d_m, k=k, H=H, NtD_modes= NTD_MODES)

    assert_allclose(I_block, I_serial, rtol=TOL, atol=TOL)


@pytest.mark.parametrize('edges', [(SEGMENT1_SIGMA_L, SEGMENT2_SIGMA_L),
                                   (SEGMENT1_SIGMA_R, SEGMENT2_SIGMA_R)])
def test_INuv(edges):
    edge_u, edge_v = edges
    I_block = block.I_Nuv(segment_u=edge_u, segment_v=edge_v, D_u=DIRECTIONS, D_v=DIRECTIONS, k=k, H=H, NtD_modes= NTD_MODES)
    I_serial = np.zeros_like(I_block)
    for n, d_n in enumerate(DIRECTIONS):
        for m, d_m in enumerate(DIRECTIONS):
            I_serial[m, n] = serial.I_Nuv(segment_u=edge_u, segment_v=edge_v, d_u=d_n, d_v=d_m, k=k, H=H, NtD_modes= NTD_MODES)

    assert_allclose(I_block, I_serial, rtol=TOL, atol=TOL)


@pytest.mark.parametrize('edges', [(SEGMENT1_SIGMA_L, SEGMENT2_SIGMA_L),
                                   (SEGMENT1_SIGMA_R, SEGMENT2_SIGMA_R)])
def test_INuNv(edges):
    edge_u, edge_v = edges
    I_block = block.I_NuNv(segment_u=edge_u, segment_v=edge_v, D_u=DIRECTIONS, D_v=DIRECTIONS, k=k, H=H, NtD_modes= NTD_MODES)
    I_serial = np.zeros_like(I_block)
    for n, d_n in enumerate(DIRECTIONS):
        for m, d_m in enumerate(DIRECTIONS):
            I_serial[m, n] = serial.I_NuNv(segment_u=edge_u, segment_v=edge_v, d_u=d_n, d_v=d_m, k=k, H=H, NtD_modes= NTD_MODES)

    assert_allclose(I_block, I_serial, rtol=TOL, atol=TOL)


@pytest.mark.parametrize('edges', [(SEGMENT1_SIGMA_L, SEGMENT2_SIGMA_L),
                                   (SEGMENT1_SIGMA_R, SEGMENT2_SIGMA_R)])
def test_INudv(edges):
    edge_u, edge_v = edges
    I_block = block.I_Nudv(segment_u=edge_u, segment_v=edge_v, D_u=DIRECTIONS, D_v=DIRECTIONS, k=k, H=H, NtD_modes= NTD_MODES)
    I_serial = np.zeros_like(I_block)
    for n, d_n in enumerate(DIRECTIONS):
        for m, d_m in enumerate(DIRECTIONS):
            I_serial[m, n] = serial.I_Nudv(segment_u=edge_u, segment_v=edge_v, d_u=d_n, d_v=d_m, k=k, H=H, NtD_modes= NTD_MODES)

    assert_allclose(I_block, I_serial, rtol=TOL, atol=TOL)