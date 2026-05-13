import numpy as np
from matplotlib.tri import Triangulation
from trefftz.numpy_types import float_array, int_array
from .core import CellLocator, TrefftzMesh


class MatplotlibLocator(CellLocator):
    def __init__(self, Tri: Triangulation):
        self.trifinder = Tri.get_trifinder()
    
    def find_cell(self, p: float_array) -> int_array | int:
        p_x, p_y = np.transpose(p)
        return self.trifinder(p_x, p_y)
    

def Mesh_from_Matplotlib(Tri: Triangulation) -> TrefftzMesh:
    '''Returns a Mesh from a Matptlotlib.tri.Triangulation'''
    points = np.column_stack([Tri.x, Tri.y])
    edges = Tri.edges
    triangles = Tri.triangles
    locator = MatplotlibLocator(Tri=Tri)
    return TrefftzMesh(points=points, edges=edges, triangles=triangles, locator=locator, cell_sets={}, edge2triangles=[] )

