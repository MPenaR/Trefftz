from .core import Mesh
import matplotlib.pyplot as plt 
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
import numpy as np

def explore_edges(mesh: Mesh):
    fig, ax = plt.subplots()
    lw = 1
    xmin, ymin = mesh._points.min(axis=0)
    xmax, ymax = mesh._points.max(axis=0)

    N_E = mesh.n_edges

    ax.add_collection(LineCollection(np.stack([mesh.edges["P"], mesh.edges["Q"]], axis=1), 
                                      colors='k', linewidths=lw))
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axis('equal')

    e = 0
    edge = mesh.edges[e]
    px, py = edge["P"]
    qx, qy = edge["Q"]
    ax.plot([px, qx], [py, qy], "b", linewidth=2*lw )
    if edge["boundary"]:
        triangle = mesh._triangles[edge["triangles"][0]]
        A, B, C = mesh._points[triangle[0]], mesh._points[triangle[1]], mesh._points[triangle[2]]
        ax.add_patch(Polygon(np.vstack([A,B,C]), facecolor='r'))
    else:
        triangle = mesh._triangles[edge["triangles"][0]]
        A, B, C = mesh._points[triangle[0]], mesh._points[triangle[1]], mesh._points[triangle[2]]
        ax.add_patch(Polygon(np.vstack([A,B,C]), facecolor='r'))

        triangle = mesh._triangles[edge["triangles"][1]]
        A, B, C = mesh._points[triangle[0]], mesh._points[triangle[1]], mesh._points[triangle[2]]
        ax.add_patch(Polygon(np.vstack([A,B,C]), facecolor='g'))

    plot_normals = True

    if plot_normals:
        ax.quiver(mesh.edges["M"][:, 0],
                  mesh.edges["M"][:, 1],
                  mesh.edges["N"][:, 0],
                  mesh.edges["N"][:, 1], angles='xy', scale_units='xy', scale=5)


    def update_plot(e):
        ax.lines[-1].remove()
        fig.canvas.draw_idle()

        if ax.patches:
            ax.patches[-1].remove()
        if ax.patches:
            ax.patches[-1].remove()

        edge = mesh.edges[e]
        px, py = edge["P"]
        qx, qy = edge["Q"]
        ax.plot([px, qx], [py, qy], "b", linewidth=2*lw )
        ax.set_title(f'Edge number: {e}, boundary: {edge["boundary"]}')
        if edge["boundary"]:
            triangle = mesh._triangles[edge["triangles"][0]]
            A, B, C = mesh._points[triangle[0]], mesh._points[triangle[1]], mesh._points[triangle[2]]
            ax.add_patch(Polygon(np.vstack([A,B,C]), facecolor='r'))
        else:
            triangle = mesh._triangles[edge["triangles"][0]]
            A, B, C = mesh._points[triangle[0]], mesh._points[triangle[1]], mesh._points[triangle[2]]
            ax.add_patch(Polygon(np.vstack([A,B,C]), facecolor='r'))

            triangle = mesh._triangles[edge["triangles"][1]]
            A, B, C = mesh._points[triangle[0]], mesh._points[triangle[1]], mesh._points[triangle[2]]
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

# def triangle_id_tester(mesh: Mesh):
#     fig, ax = plt.subplots()
#     xmin, xmax = np.min(mesh._points[:,0]), np.max(mesh._points[:,0]) 
#     ymin, ymax = np.min(mesh._points[:,1]), np.max(mesh._points[:,1])
#     N = 200
#     x = np.linspace(xmin,xmax,N)
#     y = np.linspace(ymin,ymax,N)
#     X, Y = np.meshgrid(x,y)
#     xy = np.column_stack([X.flatten(), Y.flatten()])
#     Z = mesh.get_cell(xy).reshape(X.shape).astype(float)
#     Z[Z==-1] = np.nan
#     pc = ax.pcolormesh(X,Y,Z)
#     lw = 2
#     ax.triplot(Triangulation(x=mesh._points[:,0], y=mesh._points[:,1], triangles=mesh._triangles),linewidth=lw, color='k')
#     for e_ID in mesh.boundary_edges_list:
#         P, Q = mesh._edges[e_ID]
#         p_x, p_y = mesh._points[P, :]
#         q_x, q_y = mesh._points[Q, :]
#         ax.plot([p_x,q_x],[p_y,q_y],'r',linewidth=lw)
#     fig.colorbar(mappable=pc)
#     ax.axis('equal')
#     def hover_text(x, y):
#         return f'ID = {mesh.get_cell([x, y])}'
#     annot = ax.annotate(
#         "",
#         xy=(0, 0),
#         xytext=(10, 10),
#         textcoords="offset points",
#         bbox=dict(boxstyle="round", fc="w"),
#         arrowprops=dict(arrowstyle="->"),
#     )
#     annot.set_visible(False)
#     annot.arrowprops = None
#     def on_move(event):
#         if event.inaxes != ax:
#             annot.set_visible(False)
#             fig.canvas.draw_idle()
#             return
#         xdata, ydata = event.xdata, event.ydata
#         # Update annotation
#         annot.xy = (xdata, ydata)
#         annot.set_text(hover_text(xdata, ydata))
#         annot.set_visible(True)
#         fig.canvas.draw_idle()

#     # Connect event
#     fig.canvas.mpl_connect("motion_notify_event", on_move)
#     plt.show()

