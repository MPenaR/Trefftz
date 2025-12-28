from .core import Mesh
import matplotlib.pyplot as plt 
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
import numpy as np

def _visually_test_edges(M: Mesh):
    fig, ax = plt.subplots()
    lw = 1
    xmin, ymin = M._points.min(axis=0)
    xmax, ymax = M._points.max(axis=0)

    N_E = M.n_edges

    ax.add_collection(LineCollection(np.stack([M.edges["P"], M.edges["Q"]], axis=1), 
                                      colors='k', linewidths=lw))
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axis('equal')

    e = 0
    edge = M.edges[e]
    px, py = edge["P"]
    qx, qy = edge["Q"]
    ax.plot([px, qx], [py, qy], "b", linewidth=2*lw )
    if edge["boundary"]:
        triangle = M._triangles[edge["triangles"][0]]
        A, B, C = M._points[triangle[0]], M._points[triangle[1]], M._points[triangle[2]]
        ax.add_patch(Polygon(np.vstack([A,B,C]), facecolor='r'))
    else:
        triangle = M._triangles[edge["triangles"][0]]
        A, B, C = M._points[triangle[0]], M._points[triangle[1]], M._points[triangle[2]]
        ax.add_patch(Polygon(np.vstack([A,B,C]), facecolor='r'))

        triangle = M._triangles[edge["triangles"][1]]
        A, B, C = M._points[triangle[0]], M._points[triangle[1]], M._points[triangle[2]]
        ax.add_patch(Polygon(np.vstack([A,B,C]), facecolor='g'))


    def update_plot(e):
        ax.lines[-1].remove()
        fig.canvas.draw_idle()

        if ax.patches:
            ax.patches[-1].remove()
        if ax.patches:
            ax.patches[-1].remove()

        edge = M.edges[e]
        px, py = edge["P"]
        qx, qy = edge["Q"]
        ax.plot([px, qx], [py, qy], "b", linewidth=2*lw )
        ax.set_title(f'Edge number: {e}, boundary: {edge["boundary"]}')
        if edge["boundary"]:
            triangle = M._triangles[edge["triangles"][0]]
            A, B, C = M._points[triangle[0]], M._points[triangle[1]], M._points[triangle[2]]
            ax.add_patch(Polygon(np.vstack([A,B,C]), facecolor='r'))
        else:
            triangle = M._triangles[edge["triangles"][0]]
            A, B, C = M._points[triangle[0]], M._points[triangle[1]], M._points[triangle[2]]
            ax.add_patch(Polygon(np.vstack([A,B,C]), facecolor='r'))

            triangle = M._triangles[edge["triangles"][1]]
            A, B, C = M._points[triangle[0]], M._points[triangle[1]], M._points[triangle[2]]
            ax.add_patch(Polygon(np.vstack([A,B,C]), facecolor='g'))


        fig.canvas.draw_idle()

    def on_key(event):
        nonlocal e
        if event.key == "up":
            e = min(N_E-1, e+1)
            update_plot(e)

        elif event.key == "down":
            e = max(0, e-1)
            update_plot(e)

        elif event.key == "escape":
            plt.close(fig)

    fig.canvas.mpl_connect("key_press_event", on_key)
    update_plot(e)


    plt.show()
