'''Module for defining a waveguide class, as it will be, I think, the most usefull'''
from trefftz.mesh import TrefftzMesh, EdgeType
from trefftz.mesh.from_pygmsh import Mesh_from_meshio
from typing import Optional, Callable, Any
from pygmsh.geo import Geometry
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from enum import IntEnum


# from scipy.sparse import csc_matrix
# from trefftz.numpy_types import float_array


class Region(IntEnum):
    GAMMA = 1
    SIGMA_L = 2
    SIGMA_R = 3
    SIGMA = 4


class Scatterer:
    pass

class Waveguide:
    def __init__(self, H: float = 1., R: float = 5., lc: float = 0.3,
                 scatterers: Optional[tuple[Scatterer]] = None):
        self.R = R
        self.H = H
        with Geometry() as geom:
            p0 = geom.add_point([-R, 0.], mesh_size=lc)
            p1 = geom.add_point([ R, 0.], mesh_size=lc)
            p2 = geom.add_point([ R, H ], mesh_size=lc)
            p3 = geom.add_point([-R, H ], mesh_size=lc)

            bottom = geom.add_line( p0, p1)
            right  = geom.add_line(p1, p2)
            top    = geom.add_line(p2, p3)
            left   = geom.add_line(p3, p0)

            boundary = geom.add_curve_loop([bottom, right, top, left])

            if scatterers:
                pass


            domain = geom.add_plane_surface(boundary)

            geom.add_physical(domain, label="Omega")
            geom.add_physical([bottom, top], label="Gamma")
            geom.add_physical(left, label="S_L")
            geom.add_physical(right, label="S_R")
            geom.add_physical([left, right], label="S")
            self._domain = Mesh_from_meshio(geom.generate_mesh())
    
    @property
    def domain(self) -> "TrefftzMesh":
        return self._domain
    

    def plot(self, figsize: Optional[tuple[int, int]] = (16, 2), line_width: Optional[int] = 1):
        from matplotlib.collections import LineCollection
        _, ax = plt.subplots(figsize=figsize)
        # ax.triplot(Triangulation(x=M._points[:,0], y=M._points[:,1], triangles=M._triangles),linewidth=lw, color='k')

        S = self.domain._cell_sets["S"]["line"]
        G = self.domain._cell_sets["Gamma"]["line"]  # it allows for multidimensional subsets

        inner = np.where(self.domain.edges["type"] == EdgeType.INNER)[0]

        ax.add_collection(LineCollection(np.stack([self.domain.edges[inner]["P"],
                                                   self.domain.edges[inner]["Q"]], axis=1),
                                                   colors='k', linewidths=line_width))

        ax.add_collection(LineCollection(np.stack([self.domain.edges[S]["P"],
                                                   self.domain.edges[S]["Q"]], axis=1),
                                                   colors='r', linewidths=line_width))

        ax.add_collection(LineCollection(np.stack([self.domain.edges[G]["P"],
                                                   self.domain.edges[G]["Q"]], axis=1),
                                                   colors='b', linewidths=line_width))


        ax.axis('equal')
        ax.axis('off')
        plt.show()
    

    def plot_field(self, u: Callable[[NDArray[Any], NDArray[Any]], NDArray[Any]],
                   N: int = 100, figsize: Optional[tuple[int, int]] = (16, 2), real_part: bool = False):
        x = np.linspace(-self.R, self.R, N)
        y = np.linspace(0., self.H, N)
        X, Y = np.meshgrid(x, y)
        Z = u(X, Y)

        if real_part:
            Z = np.real(Z)

        _, ax = plt.subplots(figsize=figsize)

        ax.pcolorfast((-self.R, self.R), (0., self.H), Z)
        plt.show()






