from cases.open_space.ngsolve_geometries import AnularDomain, Regions
from trefftz.dg.boundary_conditions import NtDBC, DirichletBC
import numpy as np
import trefftz.dg.serial_fluxes as serial
import trefftz.dg.block_fluxes as block
from trefftz.dg.assemblers import SerialNumerics, BlockNumerics, SerialAssembler, BlockAssembler
import cases.open_space.block_fluxes as os_block
import cases.open_space.serial_fluxes as os_serial
from trefftz.dg.basis import LinearlySpacedBasis
from trefftz.problem import Problem
from numpy.testing import assert_allclose


R = 1.
r = 0.6
lc = 0.6
mesh = AnularDomain(R=R, r=r, verbosity=0, lc=lc)

mesh.curve_region(Regions.SIGMA, radius=R)
mesh.curve_region(Regions.d_Omega, radius=r)

data = {"d_inc": np.array([-1., 0.])}

boundary_conditions = {Regions.SIGMA: NtDBC(truncating_radius=R, data=None), Regions.d_Omega: DirichletBC(data=data)}

NtD_modes = 15
a = 0.5
b = 0.5
d2 = 0.5

s_numerics = SerialNumerics(interior_kernel= serial.UltraWeakFlux(a=a, b=b),
                          local_boundary_kernels={DirichletBC: serial.DirichletFlux(a=a, data=data), NtDBC: os_serial.NtDLocal(d_2=d2, NtD_modes=NtD_modes, data=None)},
                          nonlocal_boundary_kernels={NtDBC: os_serial.NtDNonLocal(d_2=d2, NtD_modes=NtD_modes, data = None)})


b_numerics = BlockNumerics(interior_kernel= block.UltraWeakFlux(a=a, b=b),
                          local_boundary_kernels={DirichletBC: block.DirichletFlux(a=a, data=data), NtDBC: os_block.NtDLocal(d_2=d2, NtD_modes=NtD_modes, data=None)},
                          nonlocal_boundary_kernels={NtDBC: os_block.NtDNonLocal(d_2=d2, NtD_modes=NtD_modes, data = None)})

k = 16.0
Ntheta = 3
basis = LinearlySpacedBasis(N_elements=mesh.n_triangles, k=k, N_theta=Ntheta)

s_assembler = SerialAssembler(mesh=mesh, boundary_conditions=boundary_conditions, numerics=s_numerics, basis=basis)
b_assembler = BlockAssembler(mesh=mesh, boundary_conditions=boundary_conditions, numerics=b_numerics, basis=basis)

def test_A():
    P = Problem(mesh=mesh, wavenumber=k, basis=basis, boundary_conditions=boundary_conditions, assembler=s_assembler, direct_solver=True)
    A_serial = P.assemble_LHS()
    P = Problem(mesh=mesh, wavenumber=k, basis=basis, boundary_conditions=boundary_conditions, assembler=b_assembler, direct_solver=True)
    A_block = P.assemble_LHS()
    assert_allclose(A_block.toarray(), A_serial.toarray())


def test_b():
    P = Problem(mesh=mesh, wavenumber=k, basis=basis, boundary_conditions=boundary_conditions, assembler=s_assembler, direct_solver=True)
    b_serial = P.assemble_RHS()

    
    P = Problem(mesh=mesh, wavenumber=k, basis=basis, boundary_conditions=boundary_conditions, assembler=b_assembler, direct_solver=True)
    b_block = P.assemble_RHS()

    assert_allclose(b_block, b_serial)
