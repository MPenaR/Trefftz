'''Module for defining Trefftz-DG functions'''
from trefftz.dg.basis import TrefftzBasis
from numpy.typing import DTypeLike
from trefftz.mesh import TrefftzMesh
from typing import Any
from trefftz.numpy_types import complex_array, float_array
import numpy as np
from enum import Enum

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

    def __init__(self, mesh: TrefftzMesh[Any, Any], basis: TrefftzBasis, dtype: DTypeLike = np.complex128):
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
        DOFs = self.basis.T_ID_to_DOFs[T_IDs, :]  # shape (Npts, Ntheta) 

        # evaluate plane-wave basis at points
        D = self.basis.D
        k = self.basis.k # actually it should be a N_Triangles sized vector.

        xy = np.column_stack([x.ravel(), y.ravel()])

        values = np.sum(self.coefficients[DOFs] * np.exp(1j * k * np.dot(xy, D.T)), axis=1)
        values = np.where(T_IDs >=0, values, np.nan).reshape(x.shape)

        return values


class ElementwiseParameter[R: Enum, T]:
    def __init__(self, mesh: TrefftzMesh[Any, R], parameter_mapping: dict[R, T]):
        self._mesh = mesh
        self._parameter_mapping = parameter_mapping
        self._values = np.empty(mesh.n_triangles, dtype=np.asarray(list(parameter_mapping.values())).dtype)
        for reg in mesh.regions:
            self._values[mesh._reg_indexes[reg]] = parameter_mapping[reg]

    def at(self, ID: int) -> T:
        
        return np.where(ID.ravel() < 0, np.nan, self._values[ID.ravel()]).reshape(ID.shape)

    def _from_coordinates(self, x: float, y: float) -> T:
        ID = self._mesh.get_cell(x, y)
        p = self.at(ID)
        return p
    
    def __call__(self, x: float, y: float) -> T:
        return self._from_coordinates(x, y)
    
