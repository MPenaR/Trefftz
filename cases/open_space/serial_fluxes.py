from trefftz.numpy_types import float_array
from cases.open_space.serial_kernels import I_Nuv, I_Nudv, I_NuNv, I_uNv, I_Nuincv, I_Nuincdv, I_NuincNv, I_uincNv

from trefftz.dg.kernels.serial import I_uv, I_duv, I_pw_v, I_pw_dv

JACOBI_ANGER_MODES = 30
class NtDLocal:
    def __init__(self, d_2: float, NtD_modes: int, data = None):
        self.data = data
        self.d_2 = d_2
        self.NtD_modes = NtD_modes
    
    def LHS(self, edge, d_u: float_array, d_v: float_array, k: float) -> complex:
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

        I = -I_duv(edge, d_u, d_v, k) -1j*k*d_2*I_uv(edge, d_u, d_v, k)

        return I

    def RHS(self, edge, d_v: float_array, k: float) -> complex:
        plane_wave = self.data["d_inc"]
        d_inc = plane_wave.d
        d_2 = self.d_2
        NtD_modes = self.NtD_modes
        I = -I_pw_dv(edge, plane_wave, d_v, k) +  I_Nuincdv(edge, d_inc, d_v, k, NtD_modes, JACOBI_ANGER_MODES) -1j*k*d_2*(
             I_NuincNv(edge, d_inc, d_v, k, NtD_modes, JACOBI_ANGER_MODES)
            -I_Nuincv(edge, d_inc, d_v, k, NtD_modes, JACOBI_ANGER_MODES)
            -I_uincNv(edge, d_inc, d_v, k, NtD_modes, JACOBI_ANGER_MODES)
            +I_pw_v(edge, plane_wave, d_v, k))
        return I

class NtDNonLocal:
    def __init__(self, d_2: float, NtD_modes: int, data = None):
        self.data = data
        self.d_2 = d_2
        self.NtD_modes = NtD_modes
    
    def LHS(self, edge_u, edge_v, d_u: float_array, d_v: float_array, k: float) -> complex:
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

        I = I_Nudv(edge_u, edge_v, d_u, d_v, k, NtD_modes, JACOBI_ANGER_MODES) - d_2*1j*k*(I_NuNv(edge_u, edge_v, d_u, d_v, k, NtD_modes, JACOBI_ANGER_MODES)
                                                                       -I_Nuv(edge_u, edge_v, d_u, d_v, k, NtD_modes, JACOBI_ANGER_MODES)
                                                                       -I_uNv(edge_u, edge_v, d_u, d_v, k, NtD_modes, JACOBI_ANGER_MODES))
        return I
