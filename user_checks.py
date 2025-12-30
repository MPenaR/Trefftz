from trefftz.mesh.checking_utilities import explore_edges
from trefftz.mesh.from_pygmsh import CleanWaveGuide


# def Unbounded(r: float = 0.5, R: float = 5, lc: float = 0.2) -> Mesh:
#     from pygmsh.geo import Geometry
#     with Geometry() as geom:
#         center = (0, 0)

#         outer = geom.add_circle(center, R, make_surface=False, mesh_size=lc)
#         inner = geom.add_circle(center, r, make_surface=False, mesh_size=lc)

#         domain = geom.add_plane_surface(
#             outer.curve_loop,
#             holes=[inner.curve_loop]
#         )

#         geom.add_physical(domain, 'domain')
#         geom.add_physical(outer.curve_loop, 'S_R')
#         geom.add_physical(inner.curve_loop, 'S_r')

#         mesh = geom.generate_mesh()

#     return Mesh_from_meshio(mesh)

    # mesh = _triangulation_sample_mesh()
    # mesh = _gmsh_sample_mesh()
    # p = np.array([[0.2, 0.5],
    #               [0.5, 0.2]])
    # indexes = np.array([1, 0])

mesh = CleanWaveGuide(lc=0.5)
    # mesh = Unbounded()

    # _triangle_id_tester(mesh)
#     plot_waveguide(mesh, plot_tangents=True)
explore_edges(mesh)











