from trefftz.numpy_types import float_array    
from trefftz.dg.numerical_kernels import N_POINTS
import numpy as np
from cases.waveguide.numerical_NtD import NtD
from numpy import trapezoid, conj, dot

def I_uNv(segment_u, segment_v, d_u: float_array, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:

    x = segment_u["P"][0]
    # T = segment_u["T"]
    N = segment_u["N"]
    t = np.linspace(0., 1., N_POINTS)
    y = H*t
    J = H

    r = np.column_stack((x*np.ones_like(y), y))

    mask_u = (y>=segment_u["P"][1]) & (y<=segment_u["Q"][1])
    u = np.where(mask_u, np.exp(1j*k*dot(r, d_u)), 0.)

    mask_v = (y>=segment_v["P"][1]) & (y<=segment_v["Q"][1])
    v = np.where(mask_v, np.exp(1j*k*dot(r, d_v)), 0.)

    dv_dn = 1j*k*dot(d_v, N)*v
    Nv = NtD(y, dv_dn, k, H, NtD_modes)
    
    I = trapezoid(u*conj(Nv), t)*J
        
    return I

def I_Nuv(segment_u, segment_v, d_u: float_array, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:

    x = segment_u["P"][0]
    # T = segment_u["T"]
    N = segment_u["N"]
    t = np.linspace(0., 1., N_POINTS)
    y = H*t
    J = H

    r = np.column_stack((x*np.ones_like(y), y))

    mask_u = (y>=segment_u["P"][1]) & (y<=segment_u["Q"][1])
    u = np.where(mask_u, np.exp(1j*k*dot(r, d_u)), 0.)

    mask_v = (y>=segment_v["P"][1]) & (y<=segment_v["Q"][1])
    v = np.where(mask_v, np.exp(1j*k*dot(r, d_v)), 0.)

    du_dn = 1j*k*dot(d_u, N)*u
    Nu = NtD(y, du_dn, k, H, NtD_modes)
    
    I = trapezoid(Nu*conj(v), t)*J
        
    return I


def I_NuNv(segment_u, segment_v, d_u: float_array, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:

    x = segment_u["P"][0]
    # T = segment_u["T"]
    N = segment_u["N"]
    t = np.linspace(0., 1., N_POINTS)
    y = H*t
    J = H

    r = np.column_stack((x*np.ones_like(y), y))

    mask_u = (y>=segment_u["P"][1]) & (y<=segment_u["Q"][1])
    u = np.where(mask_u, np.exp(1j*k*dot(r, d_u)), 0.)

    mask_v = (y>=segment_v["P"][1]) & (y<=segment_v["Q"][1])
    v = np.where(mask_v, np.exp(1j*k*dot(r, d_v)), 0.)

    du_dn = 1j*k*dot(d_u, N)*u
    Nu = NtD(y, du_dn, k, H, NtD_modes)


    dv_dn = 1j*k*dot(d_v, N)*v
    Nv = NtD(y, dv_dn, k, H, NtD_modes)
    
    I = trapezoid(Nu*conj(Nv), t)*J
        
    return I


def I_Nudv(segment_u, segment_v, d_u: float_array, d_v: float_array, k: float, H: float, NtD_modes: int) -> complex:

    x = segment_u["P"][0]
    # T = segment_u["T"]
    N = segment_u["N"]
    t = np.linspace(0., 1., N_POINTS)
    y = H*t
    J = H

    r = np.column_stack((x*np.ones_like(y), y))

    mask_u = (y>=segment_u["P"][1]) & (y<=segment_u["Q"][1])
    u = np.where(mask_u, np.exp(1j*k*dot(r, d_u)), 0.)

    mask_v = (y>=segment_v["P"][1]) & (y<=segment_v["Q"][1])
    v = np.where(mask_v, np.exp(1j*k*dot(r, d_v)), 0.)

    du_dn = 1j*k*dot(d_u, N)*u
    Nu = NtD(y, du_dn, k, H, NtD_modes)


    dv_dn = 1j*k*dot(d_v, N)*v
    
    I = trapezoid(Nu*conj(dv_dn), t)*J
        
    return I
