import numpy as np
from trefftz.numpy_types import float_array

N_POINTS = int(1E6)


def create_numerical_functions(x: float_array, k: float, d_u: float_array, d_v: float_array, N: float_array, n_u: float | complex = 1., n_v: float | complex = 1.):
    u = np.exp(1j*k*np.dot(x, n_u*d_u))
    v = np.exp(1j*k*np.dot(x, n_v*d_v))
    du_dn = 1j*n_u*k*np.dot(N, d_u)*u
    dv_dn = 1j*n_v*k*np.dot(N, d_v)*v
    return u, v, du_dn, dv_dn


def I_uv_segment(segment,  k: float, d_u: float_array, d_v: float_array, n_u: float | complex = 1., n_v: float | complex = 1.) -> complex:
    P = segment["P"]
    T = segment["T"]
    l = segment["l"]
    N = segment["N"]

    t = np.linspace(0, 1, N_POINTS)
    x = P + l*np.outer(t, T)
    J = l

    u, v, _, _ = create_numerical_functions(x, k, d_u, d_v, N, n_u, n_v)
    I = np.trapezoid(u*np.conj(v), t)*J  # int u*conj(v) dl 
    return I



def I_uv_arc(arc, k: float, d_u: float_array, d_v: float_array, n_u: float | complex = 1., n_v: float | complex = 1.) -> complex:
    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]
    
    t = np.linspace(theta_1, theta_2, N_POINTS)
    u_r = np.column_stack([np.cos(t), np.sin(t)])
    x = R*u_r
    J = R

    u, v, _, _ = create_numerical_functions(x, k, d_u, d_v, u_r, n_u, n_v)
    I = np.trapezoid(u*np.conj(v), t)*J  # int u*conj(v) dl 
    return I


def I_duv_segment(segment, k: float, d_u: float_array, d_v: float_array, n_u: float | complex = 1., n_v: float | complex = 1.) -> complex:
    P = segment["P"]
    T = segment["T"]
    l = segment["l"]
    N = segment["N"]

    t = np.linspace(0, 1, N_POINTS)
    x = P + l*np.outer(t, T)
    J = l

    _, v, du_dn, _ = create_numerical_functions(x, k, d_u, d_v, N, n_u, n_v)
    I = np.trapezoid(du_dn*np.conj(v), t)*J  # int u*conj(v) dl 
    return I


def I_duv_arc(arc, k: float, d_u: float_array, d_v: float_array, n_u: float | complex = 1., n_v: float | complex = 1.) -> complex:
    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    t = np.linspace(theta_1, theta_2, N_POINTS)
    u_r = np.column_stack([np.cos(t), np.sin(t)])
    x = R*u_r
    J = R
    _, v, du_dn, _ = create_numerical_functions(x, k, d_u, d_v, u_r, n_u, n_v)

    I = np.trapezoid(du_dn*np.conj(v), t)*J  # int du/dn*conj(v) dl 
    return I


def I_udv_segment(segment, k: float, d_u: float_array, d_v: float_array, n_u: float | complex = 1., n_v: float | complex = 1.) -> complex:
    P = segment["P"]
    T = segment["T"]
    l = segment["l"]
    N = segment["N"]

    t = np.linspace(0, 1, N_POINTS)
    x = P + l*np.outer(t, T)
    J = l

    u, _, _, dv_dn = create_numerical_functions(x, k, d_u, d_v, N, n_u, n_v)
    I = np.trapezoid(u*np.conj(dv_dn), t)*J  # int u*conj(dv/dn) dl 
    return I


def I_udv_arc(arc, k: float, d_u: float_array, d_v: float_array, n_u: float | complex = 1., n_v: float | complex = 1.) -> complex:
    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    t = np.linspace(theta_1, theta_2, N_POINTS)
    u_r = np.column_stack([np.cos(t), np.sin(t)])
    x = R*u_r
    J = R
    u, _, _, dv_dn = create_numerical_functions(x, k, d_u, d_v, u_r, n_u, n_v)

    I = np.trapezoid(u*np.conj(dv_dn), t)*J  # int u*conj(dv/dn) dl 
    return I

def I_dudv_segment(segment, k: float, d_u: float_array, d_v: float_array, n_u: float | complex = 1., n_v: float | complex = 1.) -> complex:
    P = segment["P"]
    T = segment["T"]
    l = segment["l"]
    N = segment["N"]

    t = np.linspace(0, 1, N_POINTS)
    x = P + l*np.outer(t, T)
    J = l

    _, _, du_dn, dv_dn = create_numerical_functions(x, k, d_u, d_v, N, n_u, n_v)
    I = np.trapezoid(du_dn*np.conj(dv_dn), t)*J  # int u*conj(dv/dn) dl 
    return I


def I_dudv_arc(arc, k: float, d_u: float_array, d_v: float_array, n_u: float | complex = 1., n_v: float | complex = 1.) -> complex:
    theta_1 = arc["theta_1"]
    theta_2 = arc["theta_2"]
    R = arc["R"]

    t = np.linspace(theta_1, theta_2, N_POINTS)
    u_r = np.column_stack([np.cos(t), np.sin(t)])
    x = R*u_r
    J = R

    _, _, du_dn, dv_dn = create_numerical_functions(x, k, d_u, d_v, u_r, n_u, n_v)
    I = np.trapezoid(du_dn*np.conj(dv_dn), t)*J  # int u*conj(dv/dn) dl 
    return I
