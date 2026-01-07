'''Module for defining Trefftz-DG functions'''

from typing import Any, TYPE_CHECKING

if TYPE_CHECKING:
    from trefftz.mesh import TrefftzMesh

from trefftz.numpy_types import complex_array, float_array
import numpy as np

class ComplexFunction:
    def __init__(self, domain: "TrefftzMesh",  N_theta: int, k: float) -> None:
        self.domain = domain
        self.N_theta = N_theta
        self.N_DOF = N_theta*domain.n_triangles
        self.T_ID_to_DOFs = np.reshape(np.arange(self.N_DOF, dtype=np.int64), (-1, N_theta))
        self.k = k
        self.thetas = np.linspace(0, 2*np.pi, N_theta, endpoint=False)
        self.D = np.column_stack([np.cos(self.thetas), np.sin(self.thetas)])

    def set(self, coefficients: complex_array):
        assert len(coefficients) == self.N_DOF
        self.coefficients = coefficients


    def __call__(self, x: float | float_array, y: float | float_array) -> Any:
    
        x = np.asarray(x)
        y = np.asarray(y)
        XY = np.column_stack([x, y])
        T_IDs = self.domain.get_cell(XY)
        DOFs = self.T_ID_to_DOFs[T_IDs, :]  # (N, Ntheta)
        # z = np.sum(self.coefficients[DOFs] * np.exp(1j*self.k*(np.outer(x, self.D[:, 0]) + np.outer(y, self.D[:, 1]))), axis=1)
        z = np.sum(self.coefficients[DOFs] * np.exp(1j*self.k*np.dot(XY, np.transpose(self.D))), axis=1)
    
        return z


