from mesh.from_pygmsh import CleanWaveGuide
from mesh.plot import _visually_test_edges
    # mesh = _triangulation_sample_mesh()
    # mesh = _gmsh_sample_mesh()
    # p = np.array([[0.2, 0.5],
    #               [0.5, 0.2]])
    # indexes = np.array([1, 0])

mesh = CleanWaveGuide(lc=0.2)
    # mesh = Unbounded()

    # _triangle_id_tester(mesh)
#     plot_waveguide(mesh, plot_tangents=True)
_visually_test_edges(mesh)
