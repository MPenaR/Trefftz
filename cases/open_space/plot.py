'''Module for plotting open space problems'''

import numpy as np
import matplotlib.pyplot as plt
from trefftz.mesh.core2 import TrefftzMesh, edge_dtype, arc_dtype
def plot_openspace(mesh: TrefftzMesh, plot_tangents: bool = False, plot_normals: bool = False):
    from matplotlib.collections import LineCollection
    _, ax = plt.subplots(figsize=(8,8))
    lw = 1
    # ax.triplot(Triangulation(x=M._points[:,0], y=M._points[:,1], triangles=M._triangles),linewidth=lw, color='k')

    inner_edges = mesh.interior_edges

    ax.add_collection(LineCollection(np.stack([inner_edges["P"], inner_edges["Q"]], axis=1), colors='k', linewidths=lw))
    for reg in mesh.boundary_regions:
        edges = mesh.boundary_Edges[reg]
        if edges.dtype == arc_dtype:
            theta_1 = edges["theta_1"]
            theta_2 = edges["theta_2"]
            R = edges[0]["R"]
            O = edges[0]["O"]
            t = np.linspace(0, 1, 16)
            theta = theta_1[:, np.newaxis] +  np.outer(theta_2 - theta_1, t)
            x = O[0] + R*np.cos(theta)
            y = O[1] + R*np.sin(theta)
            ax.add_collection(LineCollection(np.stack((x, y), axis=-1), linewidths=lw))
        else:
            ax.add_collection(LineCollection(np.stack([edges["P"], edges["Q"]], axis=1), linewidths=lw))
        
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
