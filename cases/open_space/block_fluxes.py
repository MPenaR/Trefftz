from trefftz.numpy_types import float_array, complex_array
from cases.open_space.block_kernels import I_Nuv, I_Nudv, I_NuNv, I_uNv, I_Nuincv, I_Nuincdv, I_NuincNv, I_uincNv

from trefftz.dg.kernels.block import I_uv, I_duv, I_pw_v, I_dpw_v

JACOBI_ANGER_MODES = 30
class NtDLocal:
    def __init__(self, d_2: float, NtD_modes: int, data = None):
        self.data = data
        self.d_2 = d_2
        self.NtD_modes = NtD_modes
    
    def LHS(self, edge, D: float_array, n: float | complex, k: float) -> complex_array:
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

        I = -I_duv(edge, D_u=D,n_u=n, D_v=D, n_v=n, k=k) -1j*k*d_2*I_uv(edge, D_u=D,n_u=n, D_v=D, n_v=n, k=k)

        return I

    def RHS(self, edge, D: float_array, n: float | complex,  k: float) -> complex_array:
        plane_wave = self.data
        d_inc = plane_wave.d
        d_2 = self.d_2
        NtD_modes = self.NtD_modes
        I = -I_dpw_v(edge, plane_wave=plane_wave, n_v=n, D_v=D, k=k) +  I_Nuincdv(edge, d_inc, D, k, NtD_modes, JACOBI_ANGER_MODES) -1j*k*d_2*(
             I_NuincNv(edge, d_inc, D, k, NtD_modes, JACOBI_ANGER_MODES)
            -I_Nuincv(edge, d_inc, D, k, NtD_modes, JACOBI_ANGER_MODES)
            -I_uincNv(edge, d_inc, D, k, NtD_modes, JACOBI_ANGER_MODES)
            +I_pw_v(edge, plane_wave, D, n, k))

        return I

class NtDNonLocal:
    def __init__(self, d_2: float, NtD_modes: int, data = None):
        self.data = data
        self.d_2 = d_2
        self.NtD_modes = NtD_modes
    
    def LHS(self, edge_u, edge_v, D_u: float_array, n_u: float | complex, D_v: float_array, n_v: float | complex, k: float) -> complex_array:
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
        NtD_modes = self.NtD_modes

        I = I_Nudv(edge_u, edge_v, D_u, D_v, k, NtD_modes, JACOBI_ANGER_MODES) - d_2*1j*k*(I_NuNv(edge_u, edge_v, D_u, D_v, k, NtD_modes, JACOBI_ANGER_MODES)
                                                                       -I_Nuv(edge_u, edge_v, D_u, D_v, k, NtD_modes, JACOBI_ANGER_MODES)
                                                                       -I_uNv(edge_u, edge_v, D_u, D_v, k, NtD_modes, JACOBI_ANGER_MODES))
        return I
