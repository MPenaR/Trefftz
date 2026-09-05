import numpy as np
from trefftz.mesh.core import fill_arc_geometry, fill_edge_geometry, edge_dtype, arc_dtype
# tolerance 
TOL = 1E-5

# directions
N_DIRECTIONS = 3
THETAS = np.linspace(0, np.pi/2, N_DIRECTIONS, endpoint=False)
DIRECTIONS = np.column_stack([np.cos(THETAS), np.sin(THETAS)])


# generic edges
EDGE_1 = fill_edge_geometry(P=np.array([3., 3.]), Q=np.array([1., 2.]))
EDGE_2 = fill_edge_geometry(P=np.array([1., 3.]), Q=np.array([3., 2.]))
EDGES = np.array([EDGE_1, EDGE_2], dtype=edge_dtype)

# generic arc
ARC_1 = fill_arc_geometry(theta_1=30/180*np.pi, theta_2=45/180*np.pi, R=3.)
ARC_2 = fill_arc_geometry(theta_1=60/180*np.pi, theta_2=90/180*np.pi, R=3.)
ARCS = np.array([ARC_1, ARC_2], dtype=arc_dtype)


# wavenumber
k = 8.0


if __name__ == "__main__":
    print(EDGES)