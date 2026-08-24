import cases.open_space.block_kernels as block
import cases.open_space.serial_kernels as serial
import pytest
from numpy.testing import assert_allclose
from ...data import TOL, DIRECTIONS, k
from .data import ARC1_SIGMA, ARC2_SIGMA, NTD_MODES, JACOBI_ANGER_MODES
import numpy as np

@pytest.mark.slow
def test_IuNv():
    arc_u = ARC1_SIGMA
    arc_v = ARC2_SIGMA
    I_block = block.I_uNv(arc_u=arc_u, arc_v=arc_v, D_u=DIRECTIONS, D_v=DIRECTIONS, k=k, NtD_modes= NTD_MODES, JA_modes = JACOBI_ANGER_MODES)
    I_serial = np.zeros_like(I_block)
    for m, d_m in enumerate(DIRECTIONS):
        for n, d_n in enumerate(DIRECTIONS):
            I_serial[m, n] = serial.I_uNv(arc_u=arc_u, arc_v=arc_v, d_u=d_n, d_v=d_m, k=k, NtD_modes= NTD_MODES, JA_modes = JACOBI_ANGER_MODES)
    assert_allclose(I_block, I_serial, rtol=TOL, atol=TOL)

@pytest.mark.slow
def test_INuv():
    arc_u = ARC1_SIGMA
    arc_v = ARC2_SIGMA
    I_block = block.I_Nuv(arc_u=arc_u, arc_v=arc_v, D_u=DIRECTIONS, D_v=DIRECTIONS, k=k, NtD_modes= NTD_MODES, JA_modes = JACOBI_ANGER_MODES)
    I_serial = np.zeros_like(I_block)
    for m, d_m in enumerate(DIRECTIONS):
        for n, d_n in enumerate(DIRECTIONS):
            I_serial[m, n] = serial.I_Nuv(arc_u=arc_u, arc_v=arc_v, d_u=d_n, d_v=d_m, k=k, NtD_modes= NTD_MODES, JA_modes = JACOBI_ANGER_MODES)
    assert_allclose(I_block, I_serial, rtol=TOL, atol=TOL)


@pytest.mark.slow
def test_INuNv():
    arc_u = ARC1_SIGMA
    arc_v = ARC2_SIGMA
    I_block = block.I_NuNv(arc_u=arc_u, arc_v=arc_v, D_u=DIRECTIONS, D_v=DIRECTIONS, k=k, NtD_modes= NTD_MODES, JA_modes = JACOBI_ANGER_MODES)
    I_serial = np.zeros_like(I_block)
    for m, d_m in enumerate(DIRECTIONS):
        for n, d_n in enumerate(DIRECTIONS):
            I_serial[m, n] = serial.I_NuNv(arc_u=arc_u, arc_v=arc_v, d_u=d_n, d_v=d_m, k=k, NtD_modes= NTD_MODES, JA_modes = JACOBI_ANGER_MODES)
    assert_allclose(I_block, I_serial, rtol=TOL, atol=TOL)


@pytest.mark.slow
def test_INudv():
    arc_u = ARC1_SIGMA
    arc_v = ARC2_SIGMA
    I_block = block.I_Nudv(arc_u=arc_u, arc_v=arc_v, D_u=DIRECTIONS, D_v=DIRECTIONS, k=k, NtD_modes= NTD_MODES, JA_modes = JACOBI_ANGER_MODES)
    I_serial = np.zeros_like(I_block)
    for m, d_m in enumerate(DIRECTIONS):
        for n, d_n in enumerate(DIRECTIONS):
            I_serial[m, n] = serial.I_Nudv(arc_u=arc_u, arc_v=arc_v, d_u=d_n, d_v=d_m, k=k, NtD_modes= NTD_MODES, JA_modes = JACOBI_ANGER_MODES)
    assert_allclose(I_block, I_serial, rtol=TOL, atol=TOL)
