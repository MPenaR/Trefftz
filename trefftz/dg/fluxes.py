'''Module for defining flux types'''
from enum import IntEnum  #, auto

class FluxType(IntEnum):
    '''Types of boundary condition'''
    Transsmision = 0  # not really a "boundary" condition
    SoundHard = 1
    SoundSoft = 2
    Radiating = 3
    Dirichlet = 4  # when homogeneous becomes Soundsoft
    Neumann = 5    # when homogeneous becomes SoundHard
