from trefftz.numpy_types import float_array
from trefftz.dg.kernels.serial import I_uv, I_duv
from cases.waveguide.serial_kernels import I_Nuv, I_uNv, I_Nudv, I_NuNv, I_mode_v, I_mode_dv
from cases.waveguide.serial_kernels import I_Nmodedv, I_Nmodev, I_modeNv, I_NmodeNv


class NtDLocal:
    def __init__(self, d_2: float, NtD_modes: int, data=None):
        self.data = data
        self.d_2 = d_2
        self.NtD_modes = NtD_modes
    
    def LHS(self, edge, d_u: float_array, d_v: float_array, k: float) -> complex:
        r"""
        Computes the non local part of the radiating flux on vertical segment:
        

        .. math::

            \int_\Sigma \left(\right)\,\mathrm{d}\ell

        
        Parameters
        ----------
        edge_u : Edge
            Edge where u is supported.
        edge_v : Edge
            Edge where v is supported.
        d_u : (float, float) array
            Propagation direction of u.
        d_v : (float, float) array
            Propagation direction of v.
        k : float
            Wave number.
        edge : Edge
            Edge parameters.

        Returns
        -------
        I : complex
            The integral.
        """

        d_2 = self.d_2
        I = -I_duv(edge, d_u, d_v, k) -1j*k*d_2*I_uv(edge, d_u, d_v, k)
        return I

    def RHS(self, edge, d_v: float_array, k: float) -> complex:
        d_2 = self.d_2
        mode = self.data
        H = mode.H 
        NtD_modes = self.NtD_modes
        I = -I_mode_dv(edge, mode, d_v, k) + I_Nmodedv(edge, mode, d_v, k, H, NtD_modes) -1j*k*d_2*( I_NmodeNv(edge, mode, d_v, k, H, NtD_modes)
                                                                                                     -I_Nmodev(edge, mode, d_v, k, H, NtD_modes)
                                                                                                     -I_modeNv(edge, mode, d_v, k, H, NtD_modes)
                                                                                                     +I_mode_v(edge, mode, d_v, k))
        return I

# specialized one, just for testing
    # def RHS(self, edge, d_v: float_array, k: float) -> complex:
    #     d_2 = self.d_2
    #     mode = self.data
    #     NtD_modes = self.NtD_modes
    #     I = -2*I_mode_dv(edge, mode, d_v, k) + 2j*k*d_2*(I_modeNv(edge, mode, d_v, k, mode.H, NtD_modes=NtD_modes) - I_mode_v(edge, mode, d_v, k))
    #     return I


class NtDNonLocal:
    def __init__(self, H: float, d_2: float, NtD_modes: int, data=None):
        self.H = H
        self.data = data
        self.d_2 = d_2
        self.NtD_modes = NtD_modes
    
    def LHS(self, edge_u, edge_v, d_u: float_array, d_v: float_array, k: float) -> complex:
        r"""
        Computes the non local part of the radiating flux on vertical segment:
        

        .. math::

            \int_\Sigma \left(\right)\,\mathrm{d}\ell

        
        Parameters
        ----------
        edge_u : Edge
            Edge where u is supported.
        edge_v : Edge
            Edge where v is supported.
        d_u : (float, float) array
            Propagation direction of u.
        d_v : (float, float) array
            Propagation direction of v.
        k : float
            Wave number.
        edge : Edge
            Edge parameters.

        Returns
        -------
        I : complex
            The integral.
        """

        d_2 = self.d_2
        H = self.H
        NtD_modes = self.NtD_modes

        I = I_Nudv(edge_u, edge_v, d_u, d_v, k, H, NtD_modes) - d_2*1j*k*(I_NuNv(edge_u, edge_v, d_u, d_v, k, H, NtD_modes)
                                                                          -I_Nuv(edge_u, edge_v, d_u, d_v, k, H, NtD_modes)
                                                                          -I_uNv(edge_u, edge_v, d_u, d_v, k, H, NtD_modes))
        return I

    # def RHS(self, edge, d_v: float_array, k: float) -> complex:
    #     mode = self.data
    #     d_2 = self.d_2
    #     H = self.H
    #     NtD_modes = self.NtD_modes

    #     I = I_Nmodedv(edge, mode, d_v, k, H, NtD_modes) - d_2*1j*k*(I_NmodeNv(edge, mode, d_v, k, H, NtD_modes)
    #                                                                 -I_Nmodev(edge, mode, d_v, k, H, NtD_modes)
    #                                                                 -I_modeNv(edge, mode, d_v, k, H, NtD_modes))
    #     return I
