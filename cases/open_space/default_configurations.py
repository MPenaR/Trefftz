'''Here I'll try to define some default configuration to ease up the testing.'''
import numpy as np

def SoundSoft_scattered_field_problem(R: float = 1., r: float = 0.2, Lc: float = 2., lc: float = 0.05,
                          curve_scatterer: bool = False, theta_inc: float = -np.pi,
                          a: float = 0.5, b: float = 0.5, d_2: float = 0.5, NtD_modes: int = 15,
                          k: float = 8.0, N_theta: int = 14, direct_solver: bool = True):

    from cases.open_space.ngsolve_geometries import AnularDomain
    from trefftz.dg.boundary_conditions import NtDBC, DirichletBC
    from trefftz.dg.exact import PlaneWave
    from trefftz.dg.block_fluxes import UltraWeakFlux, DirichletFlux
    from cases.open_space.block_fluxes import NtDLocal, NtDNonLocal
    from trefftz.dg.assemblers import BlockNumerics
    from trefftz.dg.basis import LinearlySpacedBasis
    from trefftz.dg.parameters import Elementwise
    from trefftz.dg.assemblers import BlockAssembler
    from trefftz.problem import Problem
    from cases.open_space.exact import CircularDirichlet

    mesh = AnularDomain(R=R, r=r, Lc=Lc, lc=lc)
    Boundaries = mesh.boundaries

    if curve_scatterer: 
       mesh.curve_boundaries([Boundaries.SIGMA, Boundaries.D_OMEGA], [R, r])
    else:
       mesh.curve_boundaries([Boundaries.SIGMA], [R])

    d_inc = np.array([np.cos(theta_inc), np.sin(theta_inc)])
    data = PlaneWave(d=d_inc, A=-1., k=k)

    boundary_conditions = {Boundaries.SIGMA: NtDBC(truncating_radius=R), Boundaries.D_OMEGA: DirichletBC(data=data)}
    numerics = BlockNumerics(interior_kernel=UltraWeakFlux(a=a, b=b),
                            local_boundary_kernels={DirichletBC: DirichletFlux(a=a, data=data), NtDBC: NtDLocal(d_2=d_2, NtD_modes=NtD_modes, data=None)},
                            nonlocal_boundary_kernels={NtDBC: NtDNonLocal(d_2=d_2, NtD_modes=NtD_modes, data = None)})

    Regions = mesh.regions
    refractive_index = Elementwise(mesh, {Regions.BACKGROUND: 1., Regions.OMEGA : 1.})
    basis = LinearlySpacedBasis(N_elements=mesh.n_triangles, refractive_index=refractive_index, k=k, N_theta=N_theta)
    assembler = BlockAssembler(mesh=mesh, boundary_conditions=boundary_conditions, numerics=numerics, basis=basis)

    P = Problem(mesh=mesh,
            wavenumber=k,
            basis=basis,
            boundary_conditions=boundary_conditions,
            assembler=assembler,
            direct_solver=direct_solver,
            u=lambda x, y: CircularDirichlet(x, y, R=r, k=k, theta_inc=theta_inc))
    return P


def SoundSoft_total_field_problem(R: float = 1., r: float = 0.2, Lc: float = 2., lc: float = 0.05,
                          curve_scatterer: bool = False, theta_inc: float = -np.pi,
                          a: float = 0.5, b: float = 0.5, d_2: float = 0.5, NtD_modes: int = 15,
                          k: float = 8.0, N_theta: int = 14, direct_solver: bool = True):

    from cases.open_space.ngsolve_geometries import AnularDomain
    from trefftz.dg.boundary_conditions import NtDBC, DirichletBC
    from trefftz.dg.exact import PlaneWave
    from trefftz.dg.block_fluxes import UltraWeakFlux, DirichletFlux
    from cases.open_space.block_fluxes import NtDLocal, NtDNonLocal
    from trefftz.dg.assemblers import BlockNumerics
    from trefftz.dg.basis import LinearlySpacedBasis
    from trefftz.dg.parameters import Elementwise
    from trefftz.dg.assemblers import BlockAssembler
    from trefftz.problem import Problem
    from cases.open_space.exact import CircularDirichlet

    mesh = AnularDomain(R=R, r=r, Lc=Lc, lc=lc)
    Boundaries = mesh.boundaries

    if curve_scatterer: 
       mesh.curve_boundaries([Boundaries.SIGMA, Boundaries.D_OMEGA], [R, r])
    else:
       mesh.curve_boundaries([Boundaries.SIGMA], [R])

    d_inc = np.array([np.cos(theta_inc), np.sin(theta_inc)])
    data = PlaneWave(d=d_inc, A=1., k=k)

    boundary_conditions = {Boundaries.SIGMA: NtDBC(truncating_radius=R, data=data), Boundaries.D_OMEGA: DirichletBC(data=None)}
    numerics = BlockNumerics(interior_kernel=UltraWeakFlux(a=a, b=b),
                            local_boundary_kernels={DirichletBC: DirichletFlux(a=a, data=None), NtDBC: NtDLocal(d_2=d_2, NtD_modes=NtD_modes, data=data)},
                            nonlocal_boundary_kernels={NtDBC: NtDNonLocal(d_2=d_2, NtD_modes=NtD_modes, data = data)})

    Regions = mesh.regions
    refractive_index = Elementwise(mesh, {Regions.BACKGROUND: 1., Regions.OMEGA : 1.})
    basis = LinearlySpacedBasis(N_elements=mesh.n_triangles, refractive_index=refractive_index, k=k, N_theta=N_theta)
    assembler = BlockAssembler(mesh=mesh, boundary_conditions=boundary_conditions, numerics=numerics, basis=basis)

    def u_tot(x, y):
      u_s = CircularDirichlet(x, y, R=r, k=k, theta_inc=theta_inc)
      u_i = np.exp(1j*k*(d_inc[0]*x + d_inc[1]*y))
      RHO = np.sqrt(x**2 + y**2)
      return np.where(RHO<r, np.nan, u_s+u_i)

    P = Problem(mesh=mesh,
            wavenumber=k,
            basis=basis,
            boundary_conditions=boundary_conditions,
            assembler=assembler,
            direct_solver=direct_solver,
            u=lambda x, y: u_tot(x, y))
    return P


def Penetrable_total_field_problem(R: float = 1., r: float = 0.2, Lc: float = 2., lc: float = 0.05,
                          eps_r: float = 4., curve_scatterer: bool = False, theta_inc: float = -np.pi,
                          a: float = 0.5, b: float = 0.5, d_2: float = 0.5, NtD_modes: int = 15,
                          k: float = 8.0, N_theta: int = 14, direct_solver: bool = True):

    from cases.open_space.ngsolve_geometries import AnularDomain
    from trefftz.dg.boundary_conditions import NtDBC
    from trefftz.dg.exact import PlaneWave
    from trefftz.dg.block_fluxes import UltraWeakFlux
    from cases.open_space.block_fluxes import NtDLocal, NtDNonLocal
    from trefftz.dg.assemblers import BlockNumerics
    from trefftz.dg.basis import LinearlySpacedBasis
    from trefftz.dg.parameters import Elementwise
    from trefftz.dg.assemblers import BlockAssembler
    from trefftz.problem import Problem
    from cases.open_space.exact import nf_diel_cylinder_plane_wave

    mesh = AnularDomain(R=R, r=r, Lc=Lc, lc=lc)
    Boundaries = mesh.boundaries

    if curve_scatterer: 
       mesh.curve_boundaries([Boundaries.SIGMA, Boundaries.D_OMEGA], [R, r])
    else:
       mesh.curve_boundaries([Boundaries.SIGMA], [R])

    d_inc = np.array([np.cos(theta_inc), np.sin(theta_inc)])
    data = PlaneWave(d=d_inc, A=1., k=k)

    boundary_conditions = {Boundaries.SIGMA: NtDBC(truncating_radius=R, data=data)}
    numerics = BlockNumerics(interior_kernel=UltraWeakFlux(a=a, b=b),
                            local_boundary_kernels={NtDBC: NtDLocal(d_2=d_2, NtD_modes=NtD_modes, data=data)},
                            nonlocal_boundary_kernels={NtDBC: NtDNonLocal(d_2=d_2, NtD_modes=NtD_modes, data = data)})

    Regions = mesh.regions
    refractive_index = Elementwise(mesh, {Regions.BACKGROUND: 1., Regions.OMEGA : np.sqrt(eps_r)})
    basis = LinearlySpacedBasis(N_elements=mesh.n_triangles, refractive_index=refractive_index, k=k, N_theta=N_theta)
    assembler = BlockAssembler(mesh=mesh, boundary_conditions=boundary_conditions, numerics=numerics, basis=basis)

    def u_tot(x, y):
      u_s = nf_diel_cylinder_plane_wave(x, y, k, theta_inc=theta_inc, n=np.sqrt(eps_r), c=np.array([0., 0.]), R=r)
      u_i = np.exp(1j*k*(d_inc[0]*x + d_inc[1]*y))
      RHO = np.sqrt(x**2 + y**2)
      return np.where(RHO<r, u_s, u_s+u_i)

    P = Problem(mesh=mesh,
            wavenumber=k,
            basis=basis,
            boundary_conditions=boundary_conditions,
            assembler=assembler,
            direct_solver=direct_solver,
            u=lambda x, y: u_tot(x, y))
    return P
