import matplotlib.pyplot as plt 
from .core import Mesh
import numpy as np


def plot_waveguide(mesh: Mesh, plot_tangents: bool = False, plot_normals: bool = False):
    from matplotlib.collections import LineCollection
    _, ax = plt.subplots()

    lw = 1
    # ax.triplot(Triangulation(x=M._points[:,0], y=M._points[:,1], triangles=M._triangles),linewidth=lw, color='k')

    S = mesh._cell_sets["S"]["line"]
    G = mesh._cell_sets["Gamma"]["line"]  # it allows for multidimensional subsets

    inner = np.where(np.logical_not(mesh.edges["boundary"]))[0]

    ax.add_collection(LineCollection(np.stack([mesh.edges[inner]["P"], mesh.edges[inner]["Q"]], axis=1), 
                                     colors='k', linewidths=lw))
    ax.add_collection(LineCollection(np.stack([mesh.edges[S]["P"], mesh.edges[S]["Q"]], axis=1), 
                                     colors='r', linewidths=lw))
    ax.add_collection(LineCollection(np.stack([mesh.edges[G]["P"], mesh.edges[G]["Q"]], axis=1), 
                                     colors='b', linewidths=lw))

    if plot_tangents:
        ax.quiver(mesh.edges["M"][:, 0],
                  mesh.edges["M"][:, 1],
                  mesh.edges["T"][:, 0],
                  mesh.edges["T"][:, 1], angles='xy', scale_units='xy', scale=5)

    if plot_normals:
        ax.quiver(mesh.edges["M"][:, 0],
                  mesh.edges["M"][:, 1],
                  mesh.edges["N"][:, 0],
                  mesh.edges["N"][:, 1], angles='xy', scale_units='xy', scale=5)

    ax.axis('equal')
    plt.show()
