'''Module for defining different boundary conditions'''
from enum import Enum, auto

class BoundaryCondition(Enum):
    '''Types of boundary condition'''
    SoundHard = auto()
    SoundSoft = auto()
    Radiating = auto()
    Dirichlet = auto()  # when homogeneous becomes Soundsoft
    Neumann = auto()    # when homogeneous becomes SoundHard

