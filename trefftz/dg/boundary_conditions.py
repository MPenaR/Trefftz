'''Module for defining different boundary conditions'''
from enum import IntEnum  #, auto
from typing import Protocol, Any
from numpy.typing import NDArray
class BoundaryConditionType(IntEnum):
    '''Types of boundary condition'''
    Inner = 0  # not really a "boundary" condition
    SoundHard = 1
    SoundSoft = 2
    Radiating = 3
    Dirichlet = 4  # when homogeneous becomes Soundsoft
    Neumann = 5    # when homogeneous becomes SoundHard


class BoundaryCondition(Protocol):
    def assemble_blocks(self) -> tuple[NDArray[Any]]:
        ...
