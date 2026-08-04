'''Module for plotting open space problems'''

import numpy as np
import matplotlib.pyplot as plt
from trefftz.mesh.core2 import TrefftzMesh
def plot_openspace(mesh: TrefftzMesh, plot_tangents: bool = False, plot_normals: bool = False):
    from matplotlib.collections import LineCollection
    _, ax = plt.subplots(figsize=(8,8))
    lw = 1
    # ax.triplot(Triangulation(x=M._points[:,0], y=M._points[:,1], triangles=M._triangles),linewidth=lw, color='k')

    inner_edges = mesh.interior_edges

    ax.add_collection(LineCollection(np.stack([inner_edges["P"], inner_edges["Q"]], axis=1), colors='k', linewidths=lw))
    # ax.add_collection(LineCollection(np.stack([mesh.edges[S]["P"], mesh.edges[S]["Q"]], axis=1),
    #                                  colors='r', linewidths=lw))
    # ax.add_collection(LineCollection(np.stack([mesh.edges[G]["P"], mesh.edges[G]["Q"]], axis=1),
    #                                  colors='b', linewidths=lw))



    # if plot_tangents:
    #     ax.quiver(mesh.edges["M"][:, 0],
    #               mesh.edges["M"][:, 1],
    #               mesh.edges["T"][:, 0],
    #               mesh.edges["T"][:, 1], angles='xy', scale_units='xy', scale=5)

    # if plot_normals:
    #     ax.quiver(mesh.edges["M"][:, 0],
    #               mesh.edges["M"][:, 1],
    #               mesh.edges["N"][:, 0],
    #               mesh.edges["N"][:, 1], angles='xy', scale_units='xy', scale=5)

    ax.axis('equal')
    ax.axis('off')
    plt.show()
