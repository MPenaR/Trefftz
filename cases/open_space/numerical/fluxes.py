class NtDLocal:
    def __init__(self, R: float, d_2: float, n: int):
        self.R = R
        self.mode_n = n
        self.d_2 = d_2

    
    def LHS(self, edge, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        d2 = self.d_2
        d_n = d_phi
        d_m = d_psi
        R = self.R

        P = edge["P"]
        Q = edge["Q"]
        theta_1 = np.atan2(P[1], P[0])
        theta_2 = np.atan2(Q[1], Q[0])

        theta = np.linspace(theta_1, theta_2, N_POINTS)
        u_r = np.column_stack([np.cos(theta), np.sin(theta)])
        N = u_r
        x = R*u_r
        u = exp(1j*k*dot(x, d_n))
        v = exp(1j*k*dot(x, d_m))
        du_dn = 1j*k*dot(N, d_n)*u
        I_uv = -1j*k*d2*Int(u*conj(v), theta)*R  # int -i*k*d_2*u*conj(v) dl 
        I_duv = -Int(du_dn*conj(v), theta)*R  # int -du/dn *conj(v) dl
        I = I_uv + I_duv
        return I

    def RHS(self) -> complex:
        raise NotImplementedError
        

class NtDNonLocal:
    def __init__(self, R: float, d_2: float, n: int, N_MODES: int):
        self.R = R
        self.mode_n = n
        self.d_2 = d_2
        self.N_MODES = N_MODES
    
    def LHS(self, arc_u, arc_v, d_phi: float_array, d_psi: float_array, k: float) -> complex:
        d2 = self.d_2
        d_n = d_phi
        d_m = d_psi
        R = self.R
        N_MODES = self.N_MODES

        P_u = arc_u["P"]
        Q_u = arc_u["Q"]
        theta_1_u = np.atan2(P_u[1], P_u[0])
        theta_2_u = np.atan2(Q_u[1], Q_u[0])
        
        P_v = arc_v["P"]
        Q_v = arc_v["Q"]
        theta_1_v = np.atan2(P_v[1], P_v[0])
        theta_2_v = np.atan2(Q_v[1], Q_v[0])

        theta = np.linspace(0, 2*np.pi, N_POINTS, endpoint=False)
        u_r = np.column_stack([np.cos(theta), np.sin(theta)])
        N = u_r
        x = R*u_r

        # u = np.where(theta_1_u < theta < theta_2_u, exp(1j*k*dot(x, d_n)), 0)
        mask_u = (theta_1_u <= theta) & (theta <= theta_2_u)
        u = np.zeros_like(theta, dtype=np.complex128)
        u[mask_u] = exp(1j*k*dot(x[mask_u, :], d_n))
        # v = np.where(theta_1_v < theta < theta_2_v, exp(1j*k*dot(x, d_m)), 0)
        mask_v = (theta_1_v <= theta) & (theta <= theta_2_v)
        v = np.zeros_like(theta, dtype=np.complex128)
        v[mask_v] = exp(1j*k*dot(x[mask_v, :], d_m))

        du_dn = 1j*k*dot(N, d_n)*u
        dv_dn = 1j*k*dot(N, d_m)*v

        NtD_du = NewmanntoDirichlet(theta, du_dn, k, R, N_MODES)
        NtD_dv = NewmanntoDirichlet(theta, dv_dn, k, R, N_MODES)
        
        I_Nudv = Int(NtD_du*conj(dv_dn), theta)*R
        I_NuNv = Int(NtD_du*conj(NtD_dv), theta)*R
        I_uNv  = Int(u*conj(NtD_dv), theta)*R
        I_Nuv  = Int(NtD_du*conj(v), theta)*R
        
        I = I_Nudv - d2*1j*k*(I_NuNv - I_Nuv - I_uNv)
        return I
