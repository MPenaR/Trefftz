import numpy as np
from trefftz.mesh.core import fill_edge_geometry

NTD_MODES = 3


# generic edge on a SIGMA boundary
R = 2.
H = 1. 

EDGE1_SIGMA_L = fill_edge_geometry(P=np.array([-R, 0.2*H]), Q=np.array([-R, 0.3*H]))
EDGE1_SIGMA_L["N"] = np.array([-1., 0])

EDGE1_SIGMA_R = fill_edge_geometry(P=np.array([ R, 0.1*H]), Q=np.array([ R, 0.2*H]))
EDGE1_SIGMA_R["N"] = np.array([ 1., 0])


EDGE2_SIGMA_L = fill_edge_geometry(P=np.array([-R, 0.1*H]), Q=np.array([-R, 0.2*H]))
EDGE2_SIGMA_L["N"] = np.array([-1., 0])

EDGE2_SIGMA_R = fill_edge_geometry(P=np.array([ R, 0.2*H]), Q=np.array([ R, 0.3*H]))
EDGE2_SIGMA_R["N"] = np.array([ 1., 0])
