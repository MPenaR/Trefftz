from numpy import dot, exp, sinc, pi

class SoundHardKernel:
    '''Serial SoundHard kernel'''
    def __init__(self, d_1: float):
        self.d_1 = d_1
    
    def LHS(self, edge, d_phi, d_psi, k: float) -> complex:
        d_1 = self.d_1
        d_m = d_psi
        d_n = d_phi

        M = edge.M
        l = edge.l
        N = edge.N
        T = edge.T

        return -1j*k*l*(1 + d_1 * dot(d_n, N))*dot(d_m, N)*exp(1j*k*dot(d_n - d_m, M)) * sinc(k*l/(2*pi)*dot(d_n-d_m, T))


class UltraWeakKernel:
    '''Transmission kernel fro the UWVF'''
    def __init__(self, a: float, b: float):
        self.a = a 
        self.b = b
    
    def LHS(self, edge, d_phi, d_psi, k: float) -> complex:
        a = self.a 
        b = self.b
        d_m = d_psi
        d_n = d_phi

        k_n = k 
        k_m = k

        M = edge.M
        N = edge.N
        T = edge.T
        l = edge.l

        # I = -1j*l/2*(2*a*k + k_n*dot(d_n, N) + k_m*dot(d_m, N) + 2*b/k*k_n*dot(d_n, N)*k_m*dot(d_m, N))*exp(1j*dot(k_n*d_n - k_m*d_m, M))*sinc(l/(2*pi)*dot(k_n*d_n - k_m*d_m,T))
        return -1j*k*l*( 1/2*( k_n/k*dot(d_n, N) + k_m/k*dot(d_m, N)) + a + b*k_n/k*k_m/k*dot(d_n, N)*dot(d_m, N))*exp(1j*dot(k_n*d_n - k_m*d_m, M))*sinc(l/(2*pi)*dot(k_n*d_n - k_m*d_m,T))
