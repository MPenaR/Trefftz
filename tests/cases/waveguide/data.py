import numpy as np
from trefftz.mesh.core import fill_edge_geometry

NTD_MODES = 3


# generic segment on a SIGMA boundary
R = 2.
H = 1. 

SEGMENT1_SIGMA_L = fill_edge_geometry(P=np.array([-R, 0.2*H]), Q=np.array([-R, 0.3*H]))
SEGMENT1_SIGMA_L["N"] = np.array([-1., 0])

SEGMENT1_SIGMA_R = fill_edge_geometry(P=np.array([ R, 0.1*H]), Q=np.array([ R, 0.2*H]))
SEGMENT1_SIGMA_R["N"] = np.array([ 1., 0])


SEGMENT2_SIGMA_L = fill_edge_geometry(P=np.array([-R, 0.1*H]), Q=np.array([-R, 0.2*H]))
SEGMENT2_SIGMA_L["N"] = np.array([-1., 0])

SEGMENT2_SIGMA_R = fill_edge_geometry(P=np.array([ R, 0.2*H]), Q=np.array([ R, 0.3*H]))
SEGMENT2_SIGMA_R["N"] = np.array([ 1., 0])
