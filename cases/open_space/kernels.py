from trefftz.numpy_types import float_array
from numpy import pi, exp, conj, atan2
from numpy.lib.scimath import sqrt
from scipy.special import jv, jvp, hankel1, h1vp

from trefftz.dg.serial_kernels import JAC_ANGER_MODES, I_uv_arc, I_duv_arc


class NtDLocal_circle:
    def __init__(self, R: float, d_2: float, n: int):
        self.R = R
        self.mode_n = n
        self.d_2 = d_2
    
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

        I = -I_duv_arc(edge, d_n, d_m, k) -1j*k*d_2*I_uv_arc(edge, d_n, d_m, k)

        return I

    def RHS(self, edge, d_psi: float_array, k: float) -> complex:
        raise NotImplementedError


def Fdudn(edge, d: float_array, k: float, t: int, N_modes: int) -> complex:
    r'''
    Fourier coefficients of du/dn
    '''
    theta_1 = edge["theta_1"]
    theta_2 = edge["theta_2"]
    R = edge["R"]

    phi = atan2(d[1], d[0])

    def primitive(theta: float) -> complex:
        I = -1j*exp(-1j*t*phi)*1j**t*(jvp(t, k*R)*theta +
                                          sum( 1j**p/(1j*p)*(jvp(p+t, k*R)*exp( 1j*p*(theta - phi))
                                                    -(-1)**t*jvp(p-t, k*R)*exp(-1j*p*(theta - phi)) ) for p in range(1, N_modes)))
        return I

    I = 1j*k*(primitive(theta=theta_2) - primitive(theta=theta_1))/sqrt(2*pi)
    return I

def Fu(edge, d: float_array, k: float, t: int, N_modes: int) -> complex:
    r'''
    Fourier coefficients of u
    '''

    theta_1 = edge["theta_1"]
    theta_2 = edge["theta_2"]
    R = edge["R"]

    phi = atan2(d[1], d[0])

    def primitive(theta: float) -> complex:
        I = 1j**t*exp(-1j*phi*t)*(jv(t, k*R)*theta + sum( 1j**p/(1j*p)*(jv(p+t, k*R)*exp( 1j*p*(theta - phi))
                                                              - (-1)**t*jv(p-t, k*R)*exp(-1j*p*(theta - phi))) for p in range(1, N_modes)))
        return I

    I = (primitive(theta=theta_2) - primitive(theta=theta_1))/sqrt(2*pi)

    return I

def NtD_coeff(edge, d: float_array, k: float, t: int, N_modes: int) -> complex:
    R = edge["R"]
    return 1/k*Fdudn(edge, d, k, t, N_modes)*hankel1(t, k*R)/h1vp(t, k*R)


def I_Nuv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, NtD_modes: int, N_modes: int) -> complex:
    R = edge_u["R"]
    I = R*sum(NtD_coeff(edge_u, d_u, k, t, N_modes)*conj(Fu(edge_v, d_v, k, t, N_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I

def I_uNv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, NtD_modes: int, N_modes: int) -> complex:
    R = edge_u["R"]
    I = R*sum(Fu(edge_u, d_u, k, t, N_modes)*conj(NtD_coeff(edge_v, d_v, k, t, N_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I

def I_NuNv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, NtD_modes: int, N_modes: int) -> complex:
    R = edge_u["R"]
    I = R*sum(NtD_coeff(edge_u, d_u, k, t, N_modes)*conj(NtD_coeff(edge_v, d_v, k, t, N_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I

def I_Nudv(edge_u, edge_v, d_u: float_array, d_v: float_array, k: float, NtD_modes: int, N_modes: int) -> complex:
    R = edge_u["R"]
    I = R*sum(NtD_coeff(edge_u, d_u, k, t, N_modes)*conj(Fdudn(edge_v, d_v, k, t, N_modes)) for t in range(-NtD_modes, NtD_modes+1))
    return I



class NtDNonLocal_circle:
    def __init__(self, R: float, d_2: float, n: int, NtD_modes: int):
        self.R = R
        self.mode_n = n
        self.d_2 = d_2
        self.NtD_modes = NtD_modes
    
    def LHS(self, edge_u, edge_v, d_phi: float_array, d_psi: float_array, k: float) -> complex:
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
        d_u = d_phi
        d_v = d_psi

        NtD_modes = self.NtD_modes

        I = I_Nudv(edge_u, edge_v, d_u, d_v, k, NtD_modes, JAC_ANGER_MODES) - d_2*1j*k*(I_NuNv(edge_u, edge_v, d_u, d_v, k, NtD_modes, JAC_ANGER_MODES)
                                                                                           -I_Nuv(edge_u, edge_v, d_u, d_v, k, NtD_modes, JAC_ANGER_MODES)
                                                                                           -I_uNv(edge_u, edge_v, d_u, d_v, k, NtD_modes, JAC_ANGER_MODES))
        return I

    def RHS(self, edge, d_psi: float_array, k: float) -> complex:
        raise NotImplementedError
