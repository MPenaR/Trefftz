import numpy as np
class TrefftzBasis:
    def __init__(self, N_elements: int, N_theta: int, k: float) -> None:
        self.N_theta = N_theta
        self.k = k
        self.N_elements = N_elements
        self._N_DOF = self.N_elements * self.N_theta

        # global numbering: triangle-major order
        self.T_ID_to_DOFs = np.arange(0, self.N_DOF, dtype=np.int64).reshape(self.N_elements, self.N_theta)

        # plane-wave directions
        thetas = np.linspace(0, 2*np.pi, N_theta, endpoint=False)
        self.D = np.column_stack([np.cos(thetas), np.sin(thetas)])

    @property
    def N_DOF(self):
        return self._N_DOF

    def dofs_on_element(self, t: int):
        return self.T_ID_to_DOFs[t, :]
        
    def global_direction(self, n: int):
        '''computes the direction of planewave corresponding global index n'''
        m = n % self.N_theta
        return self.D[m, :]
