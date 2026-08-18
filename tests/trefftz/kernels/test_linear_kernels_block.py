import pytest
from numpy.testing import assert_allclose
import numpy as np
import trefftz.dg.kernels.serial as serial
import trefftz.dg.kernels.block as block
from ...data import TOL, DIRECTIONS, EDGE_1, ARC_1,  k


def test_Iuv_segment():
    I_block = block.I_uv(edge=EDGE_1, D_u=DIRECTIONS, D_v=DIRECTIONS, k=k)
    I_serial = np.zeros_like(I_block)
    for n, d_n in enumerate(DIRECTIONS):
        for m, d_m in enumerate(DIRECTIONS):
            I_serial[m, n] = serial.I_uv(edge=EDGE_1, d_u=d_n, d_v=d_m, k=k)
    assert_allclose(I_block, I_serial, rtol=TOL, atol=TOL)

def test_Iduv_segment():
    I_block = block.I_duv(edge=EDGE_1, D_u=DIRECTIONS, D_v=DIRECTIONS, k=k)
    I_serial = np.zeros_like(I_block)
    for n, d_n in enumerate(DIRECTIONS):
        for m, d_m in enumerate(DIRECTIONS):
            I_serial[m, n] = serial.I_duv(edge=EDGE_1, d_u=d_n, d_v=d_m, k=k)
    assert_allclose(I_block, I_serial, rtol=TOL, atol=TOL)

def test_Iudv_segment():
    I_block = block.I_udv(edge=EDGE_1, D_u=DIRECTIONS, D_v=DIRECTIONS, k=k)
    I_serial = np.zeros_like(I_block)
    for n, d_n in enumerate(DIRECTIONS):
        for m, d_m in enumerate(DIRECTIONS):
            I_serial[m, n] = serial.I_udv(edge=EDGE_1, d_u=d_n, d_v=d_m, k=k)
    assert_allclose(I_block, I_serial, rtol=TOL, atol=TOL)


def test_Iuv_arc():
    I_block = block.I_uv(edge=ARC_1, D_u=DIRECTIONS, D_v=DIRECTIONS, k=k)
    I_serial = np.zeros_like(I_block)
    for n, d_n in enumerate(DIRECTIONS):
        for m, d_m in enumerate(DIRECTIONS):
            I_serial[m, n] = serial.I_uv(edge=ARC_1, d_u=d_n, d_v=d_m, k=k)
    assert_allclose(I_block, I_serial, rtol=TOL, atol=TOL)

def test_Iduv_arc():
    I_block = block.I_duv(edge=ARC_1, D_u=DIRECTIONS, D_v=DIRECTIONS, k=k)
    I_serial = np.zeros_like(I_block)
    for n, d_n in enumerate(DIRECTIONS):
        for m, d_m in enumerate(DIRECTIONS):
            I_serial[m, n] = serial.I_duv(edge=ARC_1, d_u=d_n, d_v=d_m, k=k)
    assert_allclose(I_block, I_serial, rtol=TOL, atol=TOL)

def test_Iudv_arc():
    I_block = block.I_udv(edge=ARC_1, D_u=DIRECTIONS, D_v=DIRECTIONS, k=k)
    I_serial = np.zeros_like(I_block)
    for n, d_n in enumerate(DIRECTIONS):
        for m, d_m in enumerate(DIRECTIONS):
            I_serial[m, n] = serial.I_udv(edge=ARC_1, d_u=d_n, d_v=d_m, k=k)
    assert_allclose(I_block, I_serial, rtol=TOL, atol=TOL)
