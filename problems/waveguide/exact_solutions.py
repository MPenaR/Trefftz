'''Collection of exact solutions for the waveguide'''
from typing import Callable, Any
from numpy.typing import NDArray
import numpy as np


def Mode(n: int, k: float, H: float = 1., R: float = 5. ) -> Callable[[NDArray[Any], NDArray[Any]], NDArray[Any]]:
    betaH = np.emath.sqrt((k*H)**2 - (n*np.pi)**2)
    u: Callable[[NDArray[Any], NDArray[Any]], NDArray[Any]] = lambda x, y: np.exp(1j*betaH*x/H)*np.cos(n*np.pi*y/H)
    return u
