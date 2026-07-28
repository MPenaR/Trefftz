from trefftz.mesh import TrefftzMesh
from typing import Optional
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np



class Waveguide:
    def __init__(self, mesh: TrefftzMesh, R: float = 5., H: float = 1.):
        self.R = R
        self.H = H
        self.mesh = mesh

    def plot(self, figsize: tuple[int, int] | None = None, line_width: Optional[int] = 1):

        if figsize is None:
            figsize = 2*int(2*self.R/self.H), 2

        _, ax = plt.subplots(figsize=figsize)
        regions = self.mesh.boundary_regions
        interior_edges = self.mesh.interior_edges
        ax.add_collection(LineCollection(np.stack([interior_edges["P"],
                                                   interior_edges["Q"]], axis=1),
                                                   colors='k', linewidths=line_width))

        DEFAULT_COLORS = {regions.GAMMA: 'blue',
                          regions.SIGMA_L: 'red',
                          regions.SIGMA_R: 'green'}



        for region in regions:
            edges = self.mesh.edges_on(region)
            ax.add_collection(LineCollection(np.stack([edges["P"], edges["Q"]], axis=1),
                                             colors=DEFAULT_COLORS[region], linewidths=line_width))

        title_str = ", ".join([f"{region.name}: {DEFAULT_COLORS[region]}" for region in regions])

        plot_tangents = True 
        plot_normals = True


        if plot_tangents:
            ax.quiver(self.mesh.edges["M"][:, 0],
                    self.mesh.edges["M"][:, 1],
                    self.mesh.edges["T"][:, 0],
                    self.mesh.edges["T"][:, 1], angles='xy', scale_units='xy', scale=5)

        if plot_normals:
            ax.quiver(self.mesh.edges["M"][:, 0],
                    self.mesh.edges["M"][:, 1],
                    self.mesh.edges["N"][:, 0],
                    self.mesh.edges["N"][:, 1], angles='xy', scale_units='xy', scale=5)



        ax.set_title(title_str)
        ax.axis('equal')
        ax.axis('off')
        plt.show()


def CleanWaveguide(R: float, H: float, lc: float) -> Waveguide:
    from cases.waveguide.ngsolve_geometries import Empty
    return Waveguide(Empty(H=H, R=R, lc=lc), R=R, H=H)