import numpy as np 

TOL = 1E-5
N_POINTS = int(1E5)


NTH = 3
DIRECTIONS = [(np.cos(th), np.sin(th)) for th in np.linspace(0, np.pi/2, NTH, endpoint=False)]
