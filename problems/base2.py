from trefftz.mesh import TrefftzMesh, BoundaryRegions
from trefftz.dg.basis import PlanewaveBasis
from typing import Protocol, Generic, Optional, Any
from collections.abc import Mapping
from scipy.sparse.linalg import spsolve
from scipy.sparse import coo_array, csr_array, csc_array
from dataclasses import dataclass
from trefftz.numpy_types import complex_array, float_array, int_array
import numpy as np
from trefftz.dg.functions2 import TrefftzFunction
from trefftz.dg.serial_kernels import SerialLocalKernel, SerialNonLocalKernel, SerialTransmissionKernel
from trefftz.dg.assemblers import Assembler, SerialAssembler, SerialNumerics
from trefftz.dg.boundary_conditions2 import BoundaryCondition

class ExactSolution:
    ...



class Problem(Generic[BoundaryRegions]):
    def __init__(self,
                 mesh: TrefftzMesh[BoundaryRegions],
                 wavenumber: float,
                 basis: PlanewaveBasis,
                 boundary_conditions: Mapping[BoundaryRegions, BoundaryCondition],
                 numerics: SerialNumerics,
                 assembler: Assembler,
                 u: ExactSolution | None = None ):
        
        self.mesh = mesh
        self.k = wavenumber
        self.basis = basis
        self.boundary_conditions = boundary_conditions
        self.u = u
        self.numerics = numerics
        self.assembler = SerialAssembler(self.mesh,
                                         self.boundary_conditions,
                                         self.numerics,
                                         self.basis)
        self._A: coo_array | None = None
        self._b: complex_array | None = None

        self._regions_local_kernel = [region for region, bc in self.boundary_conditions.items() if type(bc) in numerics.local_boundary_kernels]
        self._regions_nonlocal_kernel = [region for region, bc in self.boundary_conditions.items() if type(bc) in numerics.nonlocal_boundary_kernels]
        self._regions_RHS_term = [region for region, bc in self.boundary_conditions.items() if bc.data is not None]

    @property
    def regions_local_kernel(self):
        return self._regions_local_kernel

    @property
    def regions_nonlocal_kernel(self):
        return self._regions_nonlocal_kernel
    
    @property
    def regions_RHS_term(self):
        return self._regions_RHS_term

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

    @property
    def b(self) -> complex_array | None:
        return self._b


    def assemble(self):
        self.assemble_RHS()
        self.assemble_LHS()

    def solve(self) -> TrefftzFunction | None:
        DIRECT = True
        if self.assembled:
            if DIRECT:
                A = csr_array(self.A)
            u_h = TrefftzFunction(self.mesh, self.basis)
            u_h.set(coefficients=spsolve(A=A, b=self.b))
            self.u_h = u_h
            return u_h
        else:
            print("Problem not fully assembled yet.")
            return None
