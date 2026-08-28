from trefftz.dg.numerical.integrators import fekete3
from numpy.testing import assert_allclose
from netgen.geom2d import SplineGeometry
from ngsolve import Mesh
import numpy as np 
from .data import TOL

W = 3.
H = 2.
def test_fekete3():
    geo = SplineGeometry()
    geo.AddRectangle(p1=(0.0, 0.0), p2=(W, H))

    ngmesh = Mesh(geo.GenerateMesh())
    points = np.array([v.point for v in ngmesh.vertices])
    S = 0.
    for T in ngmesh.faces:
        A, B, C = [points[v.nr, :] for v in T.vertices]
        S += fekete3(lambda x, y: 1., A, B, C)
    
    assert_allclose(S, W*H, rtol=TOL, atol=TOL)