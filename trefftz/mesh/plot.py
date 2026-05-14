"""
Visualization utilities for Trefftz meshes and waveguide geometries.

This module provides plotting tools for inspecting mesh structures,
boundary partitions, and geometric orientation data associated with
Trefftz discretizations.

The visualization routines are primarily intended for:

- debugging mesh connectivity
- validating boundary tagging
- inspecting tangent and normal orientations
- visualizing waveguide geometries

Functions
---------
plot_waveguide
    Plot a waveguide mesh together with boundary subsets and optional
    edge orientation vectors.

Notes
-----
The plotting utilities rely on Matplotlib and are designed for
interactive exploratory workflows rather than publication-quality
rendering.
"""


import matplotlib.pyplot as plt 
from .core import TrefftzMesh
import numpy as np


def plot_waveguide(mesh: TrefftzMesh, plot_tangents: bool = False, plot_normals: bool = False):
    """
    Plot a waveguide mesh and its boundary partitions.

    The visualization distinguishes between:

    - inner edges
    - boundary subset ``S``
    - boundary subset ``Gamma``

    using different colors.

    Optionally, edge tangent and normal vectors can also be displayed.

    Parameters
    ----------
    mesh : TrefftzMesh
        Mesh object containing geometric and boundary information.

    plot_tangents : bool, default=False
        If ``True``, display edge tangent vectors.

    plot_normals : bool, default=False
        If ``True``, display edge normal vectors.

    Notes
    -----
    The edge coloring convention is:

    - black: inner edges
    - red: boundary subset ``S``
    - blue: boundary subset ``Gamma``

    Tangent and normal vectors are displayed using Matplotlib quiver
    plots centered at edge midpoints.
    """
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
