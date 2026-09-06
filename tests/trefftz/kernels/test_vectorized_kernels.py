from numpy.testing import assert_allclose
import numpy as np
import trefftz.dg.kernels.vectorized as vectorized
import trefftz.dg.kernels.block as block
from ...data import TOL, DIRECTIONS, EDGES, ARCS,  k
from trefftz.numpy_types import float_array

N_u: float_array = np.array([1., 1.])
N_v: float_array = np.array([1., 1.])

def test_Iuv_segment():
    I_vector =vectorized.I_uv(edges=EDGES, D_u=DIRECTIONS, n_u=N_u, D_v=DIRECTIONS, n_v=N_v, k=k) 
    I_block = np.zeros_like(I_vector)
    for e, (edge, n_u, n_v) in enumerate(zip(EDGES, N_u, N_v)):
        I_block[e, :, :] = block.I_uv(edge=edge, D_u=DIRECTIONS, n_u=n_u, D_v=DIRECTIONS, n_v=n_v, k=k)
    assert_allclose(I_vector, I_block, rtol=TOL, atol=TOL)

def test_Iduv_segment():
    I_vector =vectorized.I_duv(edges=EDGES, D_u=DIRECTIONS, n_u=N_u, D_v=DIRECTIONS, n_v=N_v, k=k) 
    I_block = np.zeros_like(I_vector)
    for e, (edge, n_u, n_v) in enumerate(zip(EDGES, N_u, N_v)):
        I_block[e, :, :] = block.I_duv(edge=edge, D_u=DIRECTIONS, n_u=n_u, D_v=DIRECTIONS, n_v=n_v, k=k)

    assert_allclose(I_vector, I_block, rtol=TOL, atol=TOL)

def test_Iudv_segment():
    I_vector =vectorized.I_udv(edges=EDGES, D_u=DIRECTIONS, n_u=N_u, D_v=DIRECTIONS, n_v=N_v, k=k) 
    I_block = np.zeros_like(I_vector)
    for e, (edge, n_u, n_v) in enumerate(zip(EDGES, N_u, N_v)):
        I_block[e, :, :] = block.I_udv(edge=edge, D_u=DIRECTIONS, n_u=n_u, D_v=DIRECTIONS, n_v=n_v, k=k)

    assert_allclose(I_vector, I_block, rtol=TOL, atol=TOL)


def test_Idudv_segment():
    I_vector =vectorized.I_dudv(edges=EDGES, D_u=DIRECTIONS, n_u=N_u, D_v=DIRECTIONS, n_v=N_v, k=k) 
    I_block = np.zeros_like(I_vector)
    for e, (edge, n_u, n_v) in enumerate(zip(EDGES, N_u, N_v)):
        I_block[e, :, :] = block.I_dudv(edge=edge, D_u=DIRECTIONS, n_u=n_u, D_v=DIRECTIONS, n_v=n_v, k=k)

    assert_allclose(I_vector, I_block, rtol=TOL, atol=TOL)


def test_Iuv_arc():
    I_vector =vectorized.I_uv(edges=ARCS, D_u=DIRECTIONS, n_u=N_u, D_v=DIRECTIONS, n_v=N_v, k=k) 
    I_block = np.zeros_like(I_vector)
    for e, (arc, n_u, n_v) in enumerate(zip(ARCS, N_u, N_v)):
        I_block[e, :, :] = block.I_uv(edge=arc, D_u=DIRECTIONS, n_u=n_u, D_v=DIRECTIONS, n_v=n_v, k=k)
    assert_allclose(I_vector, I_block, rtol=TOL, atol=TOL)

def test_Iduv_arc():
    I_vector =vectorized.I_duv(edges=ARCS, D_u=DIRECTIONS, n_u=N_u, D_v=DIRECTIONS, n_v=N_v, k=k) 
    I_block = np.zeros_like(I_vector)
    for e, (arc, n_u, n_v) in enumerate(zip(ARCS, N_u, N_v)):
        I_block[e, :, :] = block.I_duv(edge=arc, D_u=DIRECTIONS, n_u=n_u, D_v=DIRECTIONS, n_v=n_v, k=k)
    assert_allclose(I_vector, I_block, rtol=TOL, atol=TOL)

def test_Iudv_arc():
    I_vector =vectorized.I_udv(edges=ARCS, D_u=DIRECTIONS, n_u=N_u, D_v=DIRECTIONS, n_v=N_v, k=k) 
    I_block = np.zeros_like(I_vector)
    for e, (arc, n_u, n_v) in enumerate(zip(ARCS, N_u, N_v)):
        I_block[e, :, :] = block.I_udv(edge=arc, D_u=DIRECTIONS, n_u=n_u, D_v=DIRECTIONS, n_v=n_v, k=k)
    assert_allclose(I_vector, I_block, rtol=TOL, atol=TOL)
