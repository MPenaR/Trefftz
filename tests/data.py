import numpy as np
from trefftz.mesh.core2 import fill_arc_geometry, fill_edge_geometry

# tolerance 
TOL = 1E-5

# directions
N_DIRECTIONS = 3
THETAS = np.linspace(0, np.pi/2, N_DIRECTIONS, endpoint=False)
DIRECTIONS = np.column_stack([np.cos(THETAS), np.sin(THETAS)])


# generic edge
EDGE_1 = fill_edge_geometry(P=np.array([3., 3.]), Q=np.array([1., 2.]))

# generic arc
ARC_1 = fill_arc_geometry(theta_1=30/180*np.pi, theta_2=45/180*np.pi, R=3.)


# wavenumber
k = 8.0
