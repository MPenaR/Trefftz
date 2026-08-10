import numpy as np
from trefftz.numpy_types import float_array

N_POINTS = int(1E5)


def create_numerical_functions(x: float_array, k: float, d_u: float_array, d_v: float_array, N: float_array):
    u = np.exp(1j*k*np.dot(x, d_u))
    v = np.exp(1j*k*np.dot(x, d_v))
    du_dn = 1j*k*np.dot(N, d_u)*u
    dv_dn = 1j*k*np.dot(N, d_v)*v
    return u, v, du_dn, dv_dn


def I_uv_segment(segment, d_u: float_array, d_v: float_array, k: float) -> complex:
    P = segment["P"]
    T = segment["T"]
    l = segment["l"]
    N = segment["N"]

    t = np.linspace(0, 1, N_POINTS)
    x = P + l*np.outer(t, T)
    J = l

    u, v, _, _ = create_numerical_functions(x, k, d_u, d_v, N)
    I = np.trapezoid(u*np.conj(v), t)*J  # int u*conj(v) dl 
    return I



def I_uv_arc(arc, d_u: float_array, d_v: float_array, k: float) -> complex:
    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]
    
    t = np.linspace(theta_1, theta_2, N_POINTS)
    u_r = np.column_stack([np.cos(t), np.sin(t)])
    x = R*u_r
    J = R

    u, v, _, _ = create_numerical_functions(x, k, d_u, d_v, u_r)
    I = np.trapezoid(u*np.conj(v), t)*J  # int u*conj(v) dl 
    return I


def I_duv_segment(segment, d_u: float_array, d_v: float_array, k: float) -> complex:
    P = segment["P"]
    T = segment["T"]
    l = segment["l"]
    N = segment["N"]

    t = np.linspace(0, 1, N_POINTS)
    x = P + l*np.outer(t, T)
    J = l

    _, v, du_dn, _ = create_numerical_functions(x, k, d_u, d_v, N)
    I = np.trapezoid(du_dn*np.conj(v), t)*J  # int u*conj(v) dl 
    return I


def I_duv_arc(arc, d_u: float_array, d_v: float_array, k: float) -> complex:
    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    t = np.linspace(theta_1, theta_2, N_POINTS)
    u_r = np.column_stack([np.cos(t), np.sin(t)])
    x = R*u_r
    J = R
    _, v, du_dn, _ = create_numerical_functions(x, k, d_u, d_v, u_r)

    I = np.trapezoid(du_dn*np.conj(v), t)*J  # int du/dn*conj(v) dl 
    return I


def I_udv_segment(segment, d_u: float_array, d_v: float_array, k: float) -> complex:
    P = segment["P"]
    T = segment["T"]
    l = segment["l"]
    N = segment["N"]

    t = np.linspace(0, 1, N_POINTS)
    x = P + l*np.outer(t, T)
    J = l

    u, _, _, dv_dn = create_numerical_functions(x, k, d_u, d_v, N)
    I = np.trapezoid(u*np.conj(dv_dn), t)*J  # int u*conj(v) dl 
    return I


def I_udv_arc(arc, d_u: float_array, d_v: float_array, k: float) -> complex:
    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    t = np.linspace(theta_1, theta_2, N_POINTS)
    u_r = np.column_stack([np.cos(t), np.sin(t)])
    x = R*u_r
    J = R
    u, _, _, dv_dn = create_numerical_functions(x, k, d_u, d_v, u_r)

    I = np.trapezoid(u*np.conj(dv_dn), t)*J  # int du/dn*conj(v) dl 
    return I
