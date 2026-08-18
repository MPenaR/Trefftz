from trefftz.numpy_types import float_array    
from trefftz.dg.numerical_kernels import N_POINTS, create_numerical_functions
import numpy as np


def I_uNv(segment_u, segment_v, d_u: float_array, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:


    P = segment["P"]
    T = segment["T"]
    l = segment["l"]
    N = segment["N"]

    t = np.linspace(0, 1, N_POINTS)
    x = P + l*np.outer(t, T)
    J = l

    u, v, _, _ = create_numerical_functions(x, k, d_u, d_v, N)
    I = np.trapezoid(u*np.conj(v), t)*J  # int u*conj(v) dl 



        P = edge_u["P"]

        t = np.linspace(0, 1, N_POINTS)
        N = edge_u["N"]
        T = edge_u["T"]
        l = edge_u["l"]
        x = P + np.outer(t,T)*l
        
        P_v = edge["P"]

        t = np.linspace(0, 1, N_POINTS)
        N = edge["N"]
        T = edge["T"]
        l = edge["l"]
        x = P + np.outer(t,T)*l

        theta = np.linspace(0, 2*np.pi, N_POINTS, endpoint=False)
        u_r = np.column_stack([np.cos(theta), np.sin(theta)])
        N = u_r
        x = R*u_r

        # u = np.where(theta_1_u < theta < theta_2_u, exp(1j*k*dot(x, d_n)), 0)
        mask_u = (theta_1_u <= theta) & (theta <= theta_2_u)
        u = np.zeros_like(theta, dtype=np.complex128)
        u[mask_u] = exp(1j*k*dot(x[mask_u, :], d_n))
        # v = np.where(theta_1_v < theta < theta_2_v, exp(1j*k*dot(x, d_m)), 0)
        mask_v = (theta_1_v <= theta) & (theta <= theta_2_v)
        v = np.zeros_like(theta, dtype=np.complex128)
        v[mask_v] = exp(1j*k*dot(x[mask_v, :], d_m))

        du_dn = 1j*k*dot(N, d_n)*u
        dv_dn = 1j*k*dot(N, d_m)*v

        NtD_du = NewmanntoDirichlet(y, du_dn, k, R, N_MODES)
        NtD_dv = NewmanntoDirichlet(y, dv_dn, k, R, N_MODES)
        
        I_Nudv = Int(NtD_du*conj(dv_dn), theta)*R
        I_NuNv = Int(NtD_du*conj(NtD_dv), theta)*R
        I_uNv  = Int(u*conj(NtD_dv), theta)*R
        I_Nuv  = Int(NtD_du*conj(v), theta)*R
        
        I = I_Nudv - d2*1j*k*(I_NuNv - I_Nuv - I_uNv)
        return I
