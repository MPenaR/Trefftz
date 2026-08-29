from trefftz.mesh import TrefftzMesh 
from trefftz.dg.basis import PlanewaveBasis
from collections.abc import Mapping
from scipy.sparse.linalg import spsolve, bicgstab
from scipy.sparse import coo_array, csr_array, csc_array
from trefftz.numpy_types import complex_array
from trefftz.dg.functions import TrefftzFunction
from trefftz.dg.assemblers import Assembler, Numerics
from trefftz.dg.boundary_conditions import BoundaryCondition
from enum import Enum
from typing import Any, Callable
from trefftz.dg.numerical.integrators import fekete3 as Int
import numpy as np

class ExactSolution:
    ...



class Problem[B: Enum, N: Numerics]:
    def __init__(self,
                 mesh: TrefftzMesh[B, Any],
                 wavenumber: float,
                 basis: PlanewaveBasis,
                 boundary_conditions: Mapping[B, BoundaryCondition],
                 assembler: Assembler[B, N],
                 u: Callable[[float, float], float] = None, 
                #  u: ExactSolution | None = None,
                 direct_solver: bool = True):
        
        self.mesh = mesh
        self.k = wavenumber
        self.basis = basis
        self.boundary_conditions = boundary_conditions
        self.u = u
        self.assembler = assembler
        self._A: coo_array | None = None
        self._b: complex_array | None = None
        self.direct_solver = direct_solver

    @property
    def N_DOF(self):
        return self.basis.N_DOF

    @property
    def assembled(self) -> bool:
        return self.A is not None and self.b is not None

    def assemble_LHS(self):
        self._A = self.assembler.assemble_LHS()
        return self._A

    @property
    def A(self) -> coo_array | None:
        return self._A

    def assemble_RHS(self):
        self._b = self.assembler.assemble_RHS()
        return self._b

    @property
    def b(self) -> complex_array | None:
        return self._b

    def assemble(self):
        self.assemble_RHS()
        self.assemble_LHS()

    def solve(self) -> TrefftzFunction | None:
        if self.assembled:
            u_h = TrefftzFunction(self.mesh, self.basis)
            if self.direct_solver:
                A = csc_array(self.A)
                coeffs = spsolve(A=A, b=self.b)
            else: 
                A = csr_array(self.A)
                coeffs, _ = bicgstab(A=A, b=self.b)
            u_h.set(coefficients=coeffs)
            self.u_h = u_h
            return u_h
        else:
            print("Problem not fully assembled yet.")
            return None
        
    def compute_error(self):
        S = 0.
        for T in self.mesh.triangles:
            A, B, C = T["A"], T["B"], T["C"]
            I = Int(lambda x, y: np.abs(self.u_h(x, y) - self.u(x, y))**2, A, B, C)
            print(I)
            S+= I
        return np.sqrt(S)
