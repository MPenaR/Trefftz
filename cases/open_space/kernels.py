from trefftz.numpy_types import float_array
from numpy import pi, dot, exp, sinc, conj, atan2, sin, cos
from numpy.linalg import norm
from numpy.lib.scimath import sqrt
from scipy.special import jv, jvp, hankel1, h1vp

# class NtDLocal_Polygonal:
#     def __init__(self, R: float, d_2: float, n: int):
#         self.R = R
#         self.mode_n = n
#         self.d_2 = d_2
#         raise NotImplementedError("Need to take care of the segment not being vertical")
    
#     def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
#         r"""
#         Computes the flux on a radiating boundary with respect to the degrees
#         of freedom from the same cell, that is:

#         TODO: it is assuming that the radiating boundary consists of a vertical segment. This should be easy to generalize.
        
#         .. math::
        
#             -\int_{E}\left(d_{2}ik\Phi_n(\mathbf{x})+\nabla \Phi_n(\mathbf{x})\cdot\mathbf{n}\right)\overline{\Psi_n(\mathbf{x})}\,\mathrm{d}S_{\mathbf{x}}

#         which can be computed as:

#         .. math::
        
#             \boxed{-ikl\left(d_{2}+\mathbf{d}_{n}\cdot\mathbf{n}\right)e^{ik\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\cdot\mathbf{M}}\mathrm{sinc}\left(\frac{kl}{2\pi}\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\mathbf{j}\right)}

#         Parameters
#         ----------
#         phi : Function
#             Trial function.
#         psi : Function
#             Test function.
#         k : float
#             Wave number.
#         edge : Edge
#             Edge parameters.
#         d_2 : float
#             Stabilyzing parameter.

#         Returns
#         -------
#         I : complex
#             The integral.
#         """

#         d_2 = self.d_2
#         d_n = d_phi
#         d_m = d_psi


#         M = edge.M
#         l = edge.l
#         N = edge.N
#         T = edge.T

#         return -1j*k*l*(d_2 + dot(d_n, N))*exp(1j*k*dot(d_n - d_m, M))*sinc(k*l/(2*pi)*dot(d_n-d_m, T))

#     def RHS(self, edge, d_psi: float_array, k: float) -> complex:

#         t = self.mode_n
#         d = d_psi
#         d_2 = self.d_2
        
#         M = edge.M
#         l = edge.l
#         N = edge.N
#         T = edge.T

#         M_x, _ = M  # I dont like it, it still assumes horizontal waveguide

#         beta = sqrt(k**2 - (t*pi/H)**2)

#         if t == 0:
#             I = 2*1j*k*l*(dot(d, N) - d_2*(1-dot(d, N)))*exp(1j*k*(M_x - dot(d, M)))*sinc(k*l/(2*pi)*dot(d, T))
#         else:
#             I = 1j*k*l*(dot(d, N) - d_2*(1-k/conj(beta)*dot(d, N)))*exp(1j*(beta*M_x - k*dot(d, M)))*( exp( 1j*t*pi*dot(M, T)/H)*sinc(k*l/(2*pi)*dot(d, T) - t*(l/(2*H)))+
#                                                                                                        exp(-1j*t*pi*dot(M, T)/H)*sinc(k*l/(2*pi)*dot(d, T) + t*(l/(2*H))))
#         return I




# class NtD_nonlocal_Polygonal:
#     def __init__(self, R: float, d_2: float, M: int): 
#         self.R = R
#         self.d_2 = d_2 
#         self.M = M
#         raise NotImplementedError("Need to take care of the segment not being vertical")

#     def LHS(self, edge_1, edge_2, d_phi: float_array, d_psi: float_array, k : float) -> complex:
#         r"""
#         Computes the flux on a radiating boundary with respect to the degrees
#         of freedom from another cell, that is:

#         TODO: it is assuming that the radiating boundary consists of a vertical segment. This should be easy to generalize.
        
#         Parameters
#         ----------
#         phi : Function
#             Trial function.
#         psi : Function
#             Test function.
#         k : float
#             Wave number.
#         edge_u : Edge
#             Edge of the triangle associated to the trial function.
#         edge_v : Edge
#             Edge of the triangle associated to the test function.
#         d_2 : float
#             Stabilyzing parameter.
#         NtD_modes : int
#             Number of modes for the approximation of the NtD map.
#         H : float
#             height of the waveguide. 

#         Returns
#         -------
#         I : complex
#             The integral.

        
#         """

#         d_2 = self.d_2
#         M = self.M

#         M_u = edge_1.M
#         l_u = edge_1.l

#         M_v = edge_2.M
#         l_v = edge_2.l


#         N = edge_1.N
#         # T = edge_1.T

#         d_n = d_phi
#         d_m = d_psi

#         I1 = -1j*k*H*dot(d_n, N)*dot(d_m,N)*d_2*exp(1j*k*(dot(d_n, M_u) - dot(d_m, M_v)))*l_u/H*l_v/H*(
#             sinc(k*l_u/(2*pi)*d_n[1])*sinc(k*l_v/(2*pi)*d_m[1]) + 1/2*sum( [ k**2 / abs(sqrt(complex(k**2 - (s*pi/H)**2)))**2 * (
#             exp( 1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] + s*l_u/(2*H)) + exp(-1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] - s*l_u/(2*H)) ) *(
#             exp(-1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] + s*l_v/(2*H)) + exp( 1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] - s*l_v/(2*H)) )
#             for s in range(1, M)]) )
        
#         I2 = -1j*k*H*dot(d_n, N)*(dot(d_m, N)-d_2)*exp(1j*k*(dot(d_n, M_u) - dot(d_m,M_v)))*l_u/H*l_v/H*(
#             sinc(k*l_u/(2*pi)*d_n[1])*sinc(k*l_v/(2*pi)*d_m[1]) + 1/2*sum( [ k / sqrt(complex(k**2 - (s*pi/H)**2)) * (
#             exp( 1j*s*pi*M_u[1]/H)*sinc(k*l_u/(2*pi)*d_n[1] + s*l_u/(2*H)) + exp(-1j*s*pi*M_u[1]/H)*sinc(k*l_u/(2*pi)*d_n[1] - s*l_u/(2*H)) ) *(
#             exp(-1j*s*pi*M_v[1]/H)*sinc(k*l_v/(2*pi)*d_m[1] + s*l_v/(2*H)) + exp( 1j*s*pi*M_v[1]/H)*sinc(k*l_v/(2*pi)*d_m[1] - s*l_v/(2*H)) )
#             for s in range(1, M)]) )
        
#         I3 = 1j*k*H*dot(d_m, N)*d_2*exp(1j*k*(dot(d_n, M_u) - dot(d_m, M_v)))*l_u/H*l_v/H*(
#             sinc(k*l_u/(2*pi)*d_n[1])*sinc(k*l_v/(2*pi)*d_m[1]) + 1/2*sum( [ k / conj(sqrt(complex(k**2 - (s*pi/H)**2))) * (
#             exp( 1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] + s*l_u/(2*H)) + exp(-1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] - s*l_u/(2*H)) ) *(
#             exp(-1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] + s*l_v/(2*H)) + exp( 1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] - s*l_v/(2*H)) )
#             for s in range(1, M)]) )

#         I = I1 + I2 + I3

#         return I


class NtDLocal_circle:
    def __init__(self, R: float, d_2: float, n: int, N_modes: int):
        self.R = R
        self.mode_n = n
        self.d_2 = d_2
        self.N_modes = N_modes
    
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        r"""
        Computes the flux on a circular radiating boundary with respect to the degrees
        of freedom from the same cell, that is:

        
        .. math::
        
            -\int_{E}\left(d_{2}ik\Phi_n(\mathbf{x})+\nabla \Phi_n(\mathbf{x})\cdot\mathbf{n}\right)\overline{\Psi_n(\mathbf{x})}\,\mathrm{d}S_{\mathbf{x}}

        TODO: old parameters
        Parameters
        ----------
        phi : Function
            Trial function.
        psi : Function
            Test function.
        k : float
            Wave number.
        edge : Edge
            Edge parameters.
        d_2 : float
            Stabilyzing parameter.

        Returns
        -------
        I : complex
            The integral.
        """

        d_2 = self.d_2
        d_n = d_phi
        d_m = d_psi
        R = self.R

        N_modes = self.N_modes

        I = -I_duv(edge=edge, d_u=d_phi, d_v=d_psi, k=k, R=R, N_modes=N_modes) - d_2*1j*k*I_uv(edge=edge, d_u=d_phi, d_v=d_psi, k=k, R=R, N_modes=N_modes)

        return I

    def RHS(self, edge, d_psi: float_array, k: float) -> complex:
        raise NotImplementedError


def I_uv(edge, d_u: float_array, d_v: float_array, k: float, R: float, N_modes: int) -> complex:
    r"""
    Computes the integral: 

    .. math::
        \int_E u\conj(v)\,\mathrm{d}\ell

    where u and v are plane waves:

    .. math::
        u(\mathbf{x}) = \exp(ik\mathbf{d}_u\cdot\mathbf{x})

    .. math::
        v(\mathbf{x}) = \exp(ik\mathbf{d}_v\cdot\mathbf{x})

    and E is an arc of a circle of radius R centered at the origin.


    Patameters
    ----------
    edge :
        Edge parameters.
    d_u : (2,) float array.
        Propagation direction of the u plane wave.
    d_v : (2,) float array.
        Propagation direction of the u plane wave.
    k : float
        Wave number.
    R : float.
        Radius of the circle defining the arc E.
    N_modes : int.
        Number of modes used in the series computation of I.
        
    Returns
    -------
    I : complex
        The integral.
    """

    P = edge["P"] 
    Q = edge["Q"]

    theta_2 = atan2(Q[1], Q[0])
    theta_1 = atan2(P[1], P[0])

    D = norm(d_u - d_v)
    phi = atan2((d_u - d_v)[1], (d_u - d_v)[0])


    I = R*(jv(0, k*R*D)*(theta_2 - theta_1) + 2*sum(1j**n*jv(n, k*R*D)/n*(sin(n*(theta_2 - phi)) - sin(n*(theta_1 - phi))) for n in range(1, N_modes)))

    return I

def I_duv(edge, d_u: float_array, d_v: float_array, k: float, R: float, N_modes: int) -> complex:
    r"""
    Computes the integral: 

    .. math::
        \int_E \nabla u\cdot\mathbf{n}\conj(v)\,\mathrm{d}\ell

    where u and v are plane waves:

    .. math::
        u(\mathbf{x}) = \exp(ik\mathbf{d}_u\cdot\mathbf{x})

    .. math::
        v(\mathbf{x}) = \exp(ik\mathbf{d}_v\cdot\mathbf{x})

    and E is an arc of a circle of radius R centered at the origin.


    Patameters
    ----------
    edge :
        Edge parameters.
    d_u : (2,) float array.
        Propagation direction of the u plane wave.
    d_v : (2,) float array.
        Propagation direction of the u plane wave.
    k : float
        Wave number.
    R : float.
        Radius of the circle defining the arc E.
    N_modes : int.
        Number of modes used in the series computation of I.
        
    Returns
    -------
    I : complex
        The integral.
    """


    P = edge["P"] 
    Q = edge["Q"]

    theta_2 = atan2(Q[1], Q[0])
    theta_1 = atan2(P[1], P[0])

    D = norm(d_u - d_v)
    phi = atan2((d_u - d_v)[1], (d_u - d_v)[0])

    phi_n = atan2(d_u[1], d_u[0])

    def primitive(theta: float) -> complex:
        I = -k*R*(jv(1, k*R*D)*cos(phi - phi_n)*(theta) +
                    sum( 1j**n/n*(- jv(n-1, k*R*D)*sin(n*(theta - phi) + (phi-phi_n)) + jv(n+1, k*R*D)*sin(n*(theta - phi) - (phi-phi_n))) for n in range(1, N_modes)))
        return I
    
    I = primitive(theta=theta_2) - primitive(theta=theta_1)

    return I


def Fdudn(edge, d: float_array, k: float, R: float, t: int, N_modes: int) -> complex:
    P = edge["P"] 
    Q = edge["Q"]

    theta_2 = atan2(Q[1], Q[0])
    theta_1 = atan2(P[1], P[0])

    phi = atan2(d[1], d[0])

    def primitive(theta: float) -> complex:
        I = -1j*exp(-1j*t*phi)*1j**t*(jvp(t, k*R)*theta +
                                          sum( 1j**p/(1j*p)*(jvp(p+t, k*R)*exp( 1j*p*(theta - phi))
                                                    -(-1)**t*jvp(p-t, k*R)*exp(-1j*p*(theta - phi)) ) for p in range(1, N_modes)))
        return I

    I = 1j*k*(primitive(theta=theta_2) - primitive(theta=theta_1))/sqrt(2*pi)
    return I

def Fu(edge, d: float_array, k: float, R: float, t: int, N_modes: int) -> complex:
    P = edge["P"] 
    Q = edge["Q"]

    theta_2 = atan2(Q[1], Q[0])
    theta_1 = atan2(P[1], P[0])

    phi = atan2(d[1], d[0])

    def primitive(theta: float) -> complex:
        I = 1j**t*exp(-1j*phi*t)*(jv(t, k*R)*theta + sum( 1j**p/(1j*p)*(jv(p+t, k*R)*exp( 1j*p*(theta - phi))
                                                              - (-1)**t*jv(p-t, k*R)*exp(-1j*p*(theta - phi))) for p in range(1, N_modes)))
        return I

    I = (primitive(theta=theta_2) - primitive(theta=theta_1))/sqrt(2*pi)

    return I

def NtD_coeff(edge, d: float_array, k: float, R: float, t: int, N_modes: int) -> complex:
    NtD_t = Fdudn(edge, d, k, R, t, N_modes)
    return 1/k*NtD_t*hankel1(t, k*R)/h1vp(t, k*R)


def I_Nuv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, R: float, NtD_modes: int, N_modes: int) -> complex:
    I = R*sum(NtD_coeff(edge_u, d_u, k, R, t, N_modes)*conj(Fu(edge_v, d_v, k, R, t, N_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I

def I_uNv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, R: float, NtD_modes: int, N_modes: int) -> complex:
    I = R*sum(Fu(edge_u, d_u, k, R, t, N_modes)*conj(NtD_coeff(edge_v, d_v, k, R, t, N_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I

def I_NuNv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, R: float, NtD_modes: int, N_modes: int) -> complex:
    I = R*sum(NtD_coeff(edge_u, d_u, k, R, t, N_modes)*conj(NtD_coeff(edge_v, d_v, k, R, t, N_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I



# class NtD_nonlocal_circle:
#     def __init__(self, R: float, d_2: float, M: int): 
#         self.R = R
#         self.d_2 = d_2 
#         self.M = M 

#     def LHS(self, edge_1: Edge, edge_2:Edge, d_phi: float_array, d_psi: float_array, k : float) -> complex:
#         r"""
#         Computes the flux on a radiating boundary with respect to the degrees
#         of freedom from another cell, that is:

#         TODO: it is assuming that the radiating boundary consists of a vertical segment. This should be easy to generalize.
        
#         Parameters
#         ----------
#         phi : Function
#             Trial function.
#         psi : Function
#             Test function.
#         k : float
#             Wave number.
#         edge_u : Edge
#             Edge of the triangle associated to the trial function.
#         edge_v : Edge
#             Edge of the triangle associated to the test function.
#         d_2 : float
#             Stabilyzing parameter.
#         NtD_modes : int
#             Number of modes for the approximation of the NtD map.
#         H : float
#             height of the waveguide. 

#         Returns
#         -------
#         I : complex
#             The integral.

        
#         """

#         d_2 = self.d_2
#         H = self.H
#         M = self.M

#         M_u = edge_1.M
#         l_u = edge_1.l

#         M_v = edge_2.M
#         l_v = edge_2.l


#         N = edge_1.N
#         # T = edge_1.T

#         d_n = d_phi
#         d_m = d_psi

#         I1 = -1j*k*H*dot(d_n, N)*dot(d_m,N)*d_2*exp(1j*k*(dot(d_n, M_u) - dot(d_m, M_v)))*l_u/H*l_v/H*(
#             sinc(k*l_u/(2*pi)*d_n[1])*sinc(k*l_v/(2*pi)*d_m[1]) + 1/2*sum( [ k**2 / abs(sqrt(complex(k**2 - (s*pi/H)**2)))**2 * (
#             exp( 1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] + s*l_u/(2*H)) + exp(-1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] - s*l_u/(2*H)) ) *(
#             exp(-1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] + s*l_v/(2*H)) + exp( 1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] - s*l_v/(2*H)) )
#             for s in range(1, M)]) )
        
#         I2 = -1j*k*H*dot(d_n, N)*(dot(d_m, N)-d_2)*exp(1j*k*(dot(d_n, M_u) - dot(d_m,M_v)))*l_u/H*l_v/H*(
#             sinc(k*l_u/(2*pi)*d_n[1])*sinc(k*l_v/(2*pi)*d_m[1]) + 1/2*sum( [ k / sqrt(complex(k**2 - (s*pi/H)**2)) * (
#             exp( 1j*s*pi*M_u[1]/H)*sinc(k*l_u/(2*pi)*d_n[1] + s*l_u/(2*H)) + exp(-1j*s*pi*M_u[1]/H)*sinc(k*l_u/(2*pi)*d_n[1] - s*l_u/(2*H)) ) *(
#             exp(-1j*s*pi*M_v[1]/H)*sinc(k*l_v/(2*pi)*d_m[1] + s*l_v/(2*H)) + exp( 1j*s*pi*M_v[1]/H)*sinc(k*l_v/(2*pi)*d_m[1] - s*l_v/(2*H)) )
#             for s in range(1, M)]) )
        
#         I3 = 1j*k*H*dot(d_m, N)*d_2*exp(1j*k*(dot(d_n, M_u) - dot(d_m, M_v)))*l_u/H*l_v/H*(
#             sinc(k*l_u/(2*pi)*d_n[1])*sinc(k*l_v/(2*pi)*d_m[1]) + 1/2*sum( [ k / conj(sqrt(complex(k**2 - (s*pi/H)**2))) * (
#             exp( 1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] + s*l_u/(2*H)) + exp(-1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] - s*l_u/(2*H)) ) *(
#             exp(-1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] + s*l_v/(2*H)) + exp( 1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] - s*l_v/(2*H)) )
#             for s in range(1, M)]) )

#         I = I1 + I2 + I3

#         return I
