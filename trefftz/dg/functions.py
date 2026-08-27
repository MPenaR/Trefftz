'''Module for defining Trefftz-DG functions'''
from trefftz.dg.basis import TrefftzBasis, PlanewaveBasis
from numpy.typing import DTypeLike
from trefftz.mesh import TrefftzMesh
from typing import Any
from trefftz.numpy_types import complex_array, float_array
import numpy as np

class ComplexFunction:
    def __init__(self, domain: TrefftzMesh[Any, Any],  N_theta: int, k: float) -> None:
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
        T_IDs = self.domain.get_cell(x, y)
        DOFs = self.T_ID_to_DOFs[T_IDs, :]  # (N, Ntheta)
        z = np.sum(self.coefficients[DOFs] * np.exp(1j*self.k*np.dot(XY, self.D.T)), axis=1)
    
        return z


class TrefftzFunction:
    '''Numerical function belonging to the finite dimensiona broken Trefftz space'''

    def __init__(self, mesh: TrefftzMesh[Any, Any], basis: PlanewaveBasis, dtype: DTypeLike = np.complex128):
        self.basis = basis
        self.mesh = mesh
        self.coefficients = np.zeros(basis.N_DOF, dtype=dtype)

    def set(self, coefficients):
        assert len(coefficients) == self.basis.N_DOF
        self.coefficients[:] = coefficients

    def __call__(self, x: float | float_array, y: float | float_array):
        x = np.asarray(x)
        y = np.asarray(y)

        T_IDs = self.mesh.get_cell(x, y).ravel()
        n = self.basis.refractive_index.at(T_IDs) # I can also do it evaluating at the points...
        DOFs = self.basis.T_ID_to_DOFs[T_IDs, :]  # shape (Npts, Ntheta) 
        # evaluate plane-wave basis at points
        D = self.basis.D
        k = self.basis.k 

        xy = np.column_stack([x.ravel(), y.ravel()])

        values = np.sum(self.coefficients[DOFs] * np.exp(1j * k * n[:, None]*np.dot(xy, D.T)), axis=1)
        values = np.where(T_IDs >=0, values, np.nan).reshape(x.shape)

        return values    
