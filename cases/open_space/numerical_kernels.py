from numpy import exp, dot, conj
from numpy import trapezoid
import numpy as np
from trefftz.numpy_types import float_array
from trefftz.dg.numerical_kernels import N_POINTS
from cases.open_space.numerical_NtD import NtD

def I_Nuv(arc_u, arc_v, d_u: float_array, d_v: float_array, k: float, NtD_modes: int) -> complex:
       
    theta_1_u = arc_u["theta_1"]
    theta_2_u = arc_u["theta_2"]
    
    theta_1_v = arc_v["theta_1"]
    theta_2_v = arc_v["theta_2"]

    R = arc_u["R"]

    theta = np.linspace(0, 2*np.pi, N_POINTS, endpoint=False)
    u_r = np.column_stack([np.cos(theta), np.sin(theta)])
    N = u_r
    x = R*u_r

    mask_u = (theta_1_u <= theta) & (theta <= theta_2_u)
    u = np.where(mask_u, exp(1j*k*dot(x, d_u)), 0.)
    
    mask_v = (theta_1_v <= theta) & (theta <= theta_2_v)
    v = np.where(mask_v, exp(1j*k*dot(x, d_v)), 0)

    du_dn = 1j*k*dot(N, d_u)*u

    Nu = NtD(theta, du_dn, k, R, NtD_modes)   
    I = trapezoid(Nu*conj(v), theta)*R
    return I

def I_uNv(arc_u, arc_v, d_u: float_array, d_v: float_array, k: float, NtD_modes: int) -> complex:
       
    theta_1_u = arc_u["theta_1"]
    theta_2_u = arc_u["theta_2"]
    
    theta_1_v = arc_v["theta_1"]
    theta_2_v = arc_v["theta_2"]

    R = arc_u["R"]

    theta = np.linspace(0, 2*np.pi, N_POINTS*80, endpoint=False)
    u_r = np.column_stack([np.cos(theta), np.sin(theta)])
    N = u_r
    x = R*u_r

    mask_u = (theta_1_u <= theta) & (theta <= theta_2_u)
    u = np.where(mask_u, exp(1j*k*dot(x, d_u)), 0.)
    
    mask_v = (theta_1_v <= theta) & (theta <= theta_2_v)
    v = np.where(mask_v, exp(1j*k*dot(x, d_v)), 0)

    dv_dn = 1j*k*dot(N, d_v)*v

    Nv = NtD(theta, dv_dn, k, R, NtD_modes)   
    I = trapezoid(u*conj(Nv), theta)*R

    return I

def num_I_NuNv(arc_u, arc_v, d_u: float_array, d_v: float_array, k: float, NtD_modes: int) -> complex:
       
    theta_1_u = arc_u["theta_1"]
    theta_2_u = arc_u["theta_2"]
    
    theta_1_v = arc_v["theta_1"]
    theta_2_v = arc_v["theta_2"]

    R = arc_u["R"]

    theta = np.linspace(0, 2*np.pi, N_POINTS*80, endpoint=False)
    u_r = np.column_stack([np.cos(theta), np.sin(theta)])
    N = u_r
    x = R*u_r

    mask_u = (theta_1_u <= theta) & (theta <= theta_2_u)
    u = np.where(mask_u, exp(1j*k*dot(x, d_u)), 0.)
    
    mask_v = (theta_1_v <= theta) & (theta <= theta_2_v)
    v = np.where(mask_v, exp(1j*k*dot(x, d_v)), 0)

    du_dn = 1j*k*dot(N, d_u)*u
    dv_dn = 1j*k*dot(N, d_v)*v

    Nu = NtD(theta, du_dn, k, R, NtD_modes)   
    Nv = NtD(theta, dv_dn, k, R, NtD_modes)   

    I = trapezoid(Nu*conj(Nv), theta)*R

    return I
