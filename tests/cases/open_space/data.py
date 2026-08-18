import numpy as np
from trefftz.mesh.core import fill_arc_geometry

NTD_MODES = 3
JACOBI_ANGER_MODES = 60

# generic arcs on a SIGMA boundary
R = 3.

ARC1_SIGMA = fill_arc_geometry(theta_1=np.pi*30/180, theta_2=np.pi*45/180, R=R)
ARC2_SIGMA = fill_arc_geometry(theta_1=np.pi*60/180, theta_2=np.pi*90/180, R=R)