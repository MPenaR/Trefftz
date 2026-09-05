import time

import cupy as cp
import numpy as np

from trefftz.dg.exact import PlaneWave
from trefftz.dg.kernels.block import linear_kernels as CPU_kernels
from trefftz.dg.kernels.GPU_block import linear_kernels as GPU_kernels


# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------

NTH = 128
N_REPEAT = 100

k = 20.0
n_u = 1.0
n_v = 1.0

theta = np.linspace(0, np.pi, NTH, endpoint=False)

D_u = np.stack(
    [np.cos(theta), np.sin(theta)],
    axis=-1,
)

D_v = D_u.copy()


# Use one of your actual segment objects here
segment = ...


# ----------------------------------------------------------------------
# NumPy
# ----------------------------------------------------------------------

# Warm-up
linear_kernels.I_uv(
    segment,
    D_u,
    n_u,
    D_v,
    n_v,
    k,
)

t0 = time.perf_counter()

for _ in range(N_REPEAT):
    result_cpu = linear_kernels.I_uv(
        segment,
        D_u,
        n_u,
        D_v,
        n_v,
        k,
    )

cpu_time = (time.perf_counter() - t0) / N_REPEAT


# ----------------------------------------------------------------------
# CuPy
# ----------------------------------------------------------------------

D_u_gpu = cp.asarray(D_u)
D_v_gpu = cp.asarray(D_v)

# Warm-up
linear_kernels_cupy.I_uv(
    segment,
    D_u_gpu,
    n_u,
    D_v_gpu,
    n_v,
    k,
)

cp.cuda.Stream.null.synchronize()

t0 = time.perf_counter()

for _ in range(N_REPEAT):
    result_gpu = linear_kernels_cupy.I_uv(
        segment,
        D_u_gpu,
        n_u,
        D_v_gpu,
        n_v,
        k,
    )

cp.cuda.Stream.null.synchronize()

gpu_time = (time.perf_counter() - t0) / N_REPEAT


# ----------------------------------------------------------------------
# Check correctness
# ----------------------------------------------------------------------

result_gpu_cpu = cp.asnumpy(result_gpu)

print(f"CPU: {cpu_time * 1e3:.3f} ms")
print(f"GPU: {gpu_time * 1e3:.3f} ms")
print(f"Speedup: {cpu_time / gpu_time:.2f}x")

print(
    "max error:",
    np.max(np.abs(result_cpu - result_gpu_cpu)),
)