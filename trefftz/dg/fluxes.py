'''Module for defining flux types'''
from enum import IntEnum
from serial_fluxes import Inner, SoundHard   # , auto

class FluxType(IntEnum):
    '''Types of fluxes'''
    TRANSMISSION = 0  # not really a "boundary" condition
    SOUNDHARD = 1
    SOUNDSOFT = 2
    RADIATING = 3
    DIRICHLET = 4  # when homogeneous becomes Soundsoft
    NEUMANN = 5    # when homogeneous becomes SoundHard


FluxKernels = {FluxType.TRANSMISSION: Inner,
               FluxType.SOUNDHARD: SoundHard}