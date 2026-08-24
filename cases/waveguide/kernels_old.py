from trefftz.numpy_types import float_array
from numpy import pi, dot, exp, sinc, conj
from numpy.lib.scimath import sqrt

class NtDLocal:
    def __init__(self, R: float, d_2: float, n: int, H: float):
        self.R = R
        self.mode_n = n
        self.H = H
        self.d_2 = d_2
    
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        r"""
        Computes the flux on a radiating boundary with respect to the degrees
        of freedom from the same cell, that is:

        TODO: it is assuming that the radiating boundary consists of a vertical segment. This should be easy to generalize.
        
        .. math::
        
            -\int_{E}\left(d_{2}ik\Phi_n(\mathbf{x})+\nabla \Phi_n(\mathbf{x})\cdot\mathbf{n}\right)\overline{\Psi_n(\mathbf{x})}\,\mathrm{d}S_{\mathbf{x}}

        which can be computed as:

        .. math::
        
            \boxed{-ikl\left(d_{2}+\mathbf{d}_{n}\cdot\mathbf{n}\right)e^{ik\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\cdot\mathbf{M}}\mathrm{sinc}\left(\frac{kl}{2\pi}\left(\mathbf{d}_{n}-\mathbf{d}_{m}\right)\mathbf{j}\right)}

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


        M = edge["M"]
        l = edge["l"]
        N = edge["N"]
        T = edge["T"]

        return -1j*k*l*(d_2 + dot(d_n, N))*exp(1j*k*dot(d_n - d_m, M))*sinc(k*l/(2*pi)*dot(d_n-d_m, T))

    def RHS(self, edge, d_psi: float_array, k: float) -> complex:

        t = self.mode_n
        d = d_psi
        d_2 = self.d_2
        H = self.H
        
        M = edge["M"]
        l = edge["l"]
        N = edge["N"]
        T = edge["T"]

        M_x, _ = M  # I dont like it, it still assumes horizontal waveguide

        beta = sqrt(k**2 - (t*pi/H)**2)

        if t == 0:
            I = 2*1j*k*l*(dot(d, N) - d_2*(1-dot(d, N)))*exp(1j*k*(M_x - dot(d, M)))*sinc(k*l/(2*pi)*dot(d, T))
        else:
            I = 1j*k*l*(dot(d, N) - d_2*(1-k/conj(beta)*dot(d, N)))*exp(1j*(beta*M_x - k*dot(d, M)))*( exp( 1j*t*pi*dot(M, T)/H)*sinc(k*l/(2*pi)*dot(d, T) - t*(l/(2*H)))+
                                                                                                       exp(-1j*t*pi*dot(M, T)/H)*sinc(k*l/(2*pi)*dot(d, T) + t*(l/(2*H))))
        return I




class NtD_nonlocal:
    def __init__(self, R: float, d_2: float, M: int, H: float): 
        self.R = R
        self.d_2 = d_2 
        self.M = M 
        self.H = H

    def LHS(self, edge_u, edge_v, d_phi: float_array, d_psi: float_array, k : float) -> complex:
        r"""
        Computes the flux on a radiating boundary with respect to the degrees
        of freedom from another cell, that is:

        TODO: it is assuming that the radiating boundary consists of a vertical segment. This should be easy to generalize.
        
        Parameters
        ----------
        phi : Function
            Trial function.
        psi : Function
            Test function.
        k : float
            Wave number.
        edge_u : Edge
            Edge of the triangle associated to the trial function.
        edge_v : Edge
            Edge of the triangle associated to the test function.
        d_2 : float
            Stabilyzing parameter.
        NtD_modes : int
            Number of modes for the approximation of the NtD map.
        H : float
            height of the waveguide. 

        Returns
        -------
        I : complex
            The integral.

        
        """

        d_2 = self.d_2
        H = self.H
        M = self.M

        M_u = edge_u["M"]
        l_u = edge_u["l"]

        M_v = edge_v["M"]
        l_v = edge_v["l"]


        N = edge_u["N"]
        # T = edge_1.T

        d_n = d_phi
        d_m = d_psi

        I1 = -1j*k*H*dot(d_n, N)*dot(d_m,N)*d_2*exp(1j*k*(dot(d_n, M_u) - dot(d_m, M_v)))*l_u/H*l_v/H*(
            sinc(k*l_u/(2*pi)*d_n[1])*sinc(k*l_v/(2*pi)*d_m[1]) + 1/2*sum( [ k**2 / abs(sqrt(complex(k**2 - (s*pi/H)**2)))**2 * (
            exp( 1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] + s*l_u/(2*H)) + exp(-1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] - s*l_u/(2*H)) ) *(
            exp(-1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] + s*l_v/(2*H)) + exp( 1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] - s*l_v/(2*H)) )
            for s in range(1, M)]) )
        
        I2 = -1j*k*H*dot(d_n, N)*(dot(d_m, N)-d_2)*exp(1j*k*(dot(d_n, M_u) - dot(d_m,M_v)))*l_u/H*l_v/H*(
            sinc(k*l_u/(2*pi)*d_n[1])*sinc(k*l_v/(2*pi)*d_m[1]) + 1/2*sum( [ k / sqrt(complex(k**2 - (s*pi/H)**2)) * (
            exp( 1j*s*pi*M_u[1]/H)*sinc(k*l_u/(2*pi)*d_n[1] + s*l_u/(2*H)) + exp(-1j*s*pi*M_u[1]/H)*sinc(k*l_u/(2*pi)*d_n[1] - s*l_u/(2*H)) ) *(
            exp(-1j*s*pi*M_v[1]/H)*sinc(k*l_v/(2*pi)*d_m[1] + s*l_v/(2*H)) + exp( 1j*s*pi*M_v[1]/H)*sinc(k*l_v/(2*pi)*d_m[1] - s*l_v/(2*H)) )
            for s in range(1, M)]) )
        
        I3 = 1j*k*H*dot(d_m, N)*d_2*exp(1j*k*(dot(d_n, M_u) - dot(d_m, M_v)))*l_u/H*l_v/H*(
            sinc(k*l_u/(2*pi)*d_n[1])*sinc(k*l_v/(2*pi)*d_m[1]) + 1/2*sum( [ k / conj(sqrt(complex(k**2 - (s*pi/H)**2))) * (
            exp( 1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] + s*l_u/(2*H)) + exp(-1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] - s*l_u/(2*H)) ) *(
            exp(-1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] + s*l_v/(2*H)) + exp( 1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] - s*l_v/(2*H)) )
            for s in range(1, M)]) )

        I = I1 + I2 + I3

        return I

# Check Implementation
def uNdv(edge_u, edge_v, d_phi: float_array, d_psi: float_array, k : float) -> complex:
    r"""
    Computes the integral: 

    .. math ::
        \int_E u\conj(\mathrm{NtD}(\grad v \cdot \mathbf{n}))
    TODO: it is assuming that the radiating boundary consists of a vertical segment. This should be easy to generalize.
    
    Parameters
    ----------
    phi : Function
        Trial function.
    psi : Function
        Test function.
    k : float
        Wave number.
    edge_u : Edge
        Edge of the triangle associated to the trial function.
    edge_v : Edge
        Edge of the triangle associated to the test function.
    d_2 : float
        Stabilyzing parameter.
    NtD_modes : int
        Number of modes for the approximation of the NtD map.
    H : float
        height of the waveguide. 

    Returns
    -------
    I : complex
        The integral.

    
    """

    d_2 = self.d_2
    H = self.H
    M = self.M

    M_u = edge_u["M"]
    l_u = edge_u["l"]

    M_v = edge_v["M"]
    l_v = edge_v["l"]


    N = edge_u["N"]
    # T = edge_1.T

    d_n = d_phi
    d_m = d_psi

    I = -1j*k*H*dot(d_n, N)*dot(d_m,N)*d_2*exp(1j*k*(dot(d_n, M_u) - dot(d_m, M_v)))*l_u/H*l_v/H*(
        sinc(k*l_u/(2*pi)*d_n[1])*sinc(k*l_v/(2*pi)*d_m[1]) + 1/2*sum( [ k**2 / abs(sqrt(complex(k**2 - (s*pi/H)**2)))**2 * (
        exp( 1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] + s*l_u/(2*H)) + exp(-1j*s*pi/H*M_u[1])*sinc(k*l_u/(2*pi)*d_n[1] - s*l_u/(2*H)) ) *(
        exp(-1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] + s*l_v/(2*H)) + exp( 1j*s*pi/H*M_v[1])*sinc(k*l_v/(2*pi)*d_m[1] - s*l_v/(2*H)) )
        for s in range(1, M)]) )

    return I