from trefftz.mesh import TrefftzMesh  #, BoundaryRegions
from trefftz.dg.basis import PlanewaveBasis
from collections.abc import Mapping
from scipy.sparse.linalg import spsolve, gmres, bicgstab
from scipy.sparse import coo_array, csr_array, csc_array
from trefftz.numpy_types import complex_array
from trefftz.dg.functions import TrefftzFunction
from trefftz.dg.assemblers import Assembler, SerialAssembler, SerialNumerics, Numerics
from trefftz.dg.boundary_conditions import BoundaryCondition
from enum import Enum


class ExactSolution:
    ...



class Problem[BR: Enum, N: Numerics]:
    def __init__(self,
                 mesh: TrefftzMesh[BR],
                 wavenumber: float,
                 basis: PlanewaveBasis,
                 boundary_conditions: Mapping[BR, BoundaryCondition],
                #  numerics: SerialNumerics,
                 assembler: Assembler[BR, N],
                 u: ExactSolution | None = None,
                 direct_solver: bool = True):
        
        self.mesh = mesh
        self.k = wavenumber
        self.basis = basis
        self.boundary_conditions = boundary_conditions
        self.u = u
        # self.numerics = numerics
        self.assembler = assembler
        self._A: coo_array | None = None
        self._b: complex_array | None = None
        self.direct_solver = direct_solver
        # self._regions_local_kernel = [region for region, bc in self.boundary_conditions.items() if type(bc) in numerics.local_boundary_kernels]
        # self._regions_nonlocal_kernel = [region for region, bc in self.boundary_conditions.items() if type(bc) in numerics.nonlocal_boundary_kernels]
        # self._regions_RHS_term = [region for region, bc in self.boundary_conditions.items() if bc.data is not None]

    # @property
    # def regions_local_kernel(self):
    #     return self._regions_local_kernel

    # @property
    # def regions_nonlocal_kernel(self):
    #     return self._regions_nonlocal_kernel
    
    # @property
    # def regions_RHS_term(self):
    #     return self._regions_RHS_term

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
        # some code like || u-u_h ||²
        raise NotImplementedError
