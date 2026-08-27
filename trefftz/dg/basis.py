"""
Plane-wave Trefftz basis functions and direction-generation utilities.

This module defines data structures and helper functions for constructing
Trefftz approximation spaces based on plane waves.

The implementation provides:

- storage and management of plane-wave directions
- degree-of-freedom indexing for element-wise bases
- utilities for generating directional distributions
- visualization tools for basis directions

The basis functions are intended for Trefftz and discontinuous Galerkin
methods applied to wave propagation problems such as the Helmholtz
equation.

Classes
-------
PlanewaveBasis
    Plane-wave basis representation with explicit directional storage.

TrefftzBasis
    Legacy Trefftz basis implementation using uniformly distributed
    plane-wave directions.

Functions
---------
LinearlySpacedBasis
    Construct a basis with uniformly distributed directions.

TayloredBasis
    Construct a waveguide-adapted basis using propagating modal
    directions.

Notes
-----
Plane-wave directions are represented as unit vectors:


::contentReference[oaicite:0]{index=0}


The global degree-of-freedom numbering follows an element-major layout,
where all basis functions associated with one element are stored
contiguously.
"""

import numpy as np
from trefftz.numpy_types import float_array, int_array
from dataclasses import dataclass
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from trefftz.dg.parameters import Elementwise
from typing import Any

@dataclass(slots=True)
class PlanewaveBasis:
    """
    Plane-wave Trefftz basis for element-wise discretizations.

    This class stores the directional information and global
    degree-of-freedom numbering associated with a plane-wave Trefftz
    approximation space.

    Each mesh element is assigned the same set of plane-wave directions,
    generating a local basis of oscillatory functions.

    Parameters
    ----------
    N_elements : int
        Number of mesh elements.

    thetas : float_array
        Array of propagation angles in radians.

    k : float
        Wavenumber.

    Attributes
    ----------
    k : float
        Wavenumber of the basis functions.

    D : float_array
        Array of unit propagation directions with shape
        ``(N_theta, 2)``.

    T_ID_to_DOFs : int_array
        Mapping from element indices to global degrees of freedom.

    Notes
    -----
    Plane-wave directions are computed as:

    :contentReference[oaicite:1]{index=1}

    The global numbering uses an element-major ordering:

    :contentReference[oaicite:2]{index=2}
    """

    k: float #not really sure if k should be an attribute. 
    refractive_index: Elementwise[Any, float]
    _N_theta: int
    _N_DOF: int
    _N_elements: int
    _T_ID_to_DOFs: int_array
    _D: float_array
    _D_D: float_array
        
    def __init__(self, N_elements: int, refractive_index: Elementwise[Any, float], thetas: float_array, k: float) -> None:
        self._N_theta = len(thetas)
        self.k = k
        self._N_elements = N_elements
        self._N_DOF = self.N_elements * self.N_theta

        # global numbering: triangle-major order
        self._T_ID_to_DOFs = np.arange(0, self.N_DOF, dtype=np.int64).reshape(self.N_elements, self.N_theta)

        # plane-wave directions
        self._D = np.column_stack([np.cos(thetas), np.sin(thetas)])
        self._D_D = self._D[:, None, :] - self._D[None, :, :]
        self.refractive_index = refractive_index

    @property
    def N_theta(self):
        "Number of directions per element"
        return self._N_theta

    @property
    def N_elements(self):
        "Number of elements"
        return self._N_elements

    @property
    def N_DOF(self):
        '''Number of degrees of freedom'''
        return self._N_DOF
    
    @property
    def T_ID_to_DOFs(self):
        '''degrees of freedon per elements (array, Nelements x N theta)'''
        return self._T_ID_to_DOFs
    
    @property
    def D(self):
        '''Ntheta x 2, array of directions D = (d_ij) with:
          d_i1 = cos(theta_i)
          d_i2 = sin(theta_i) '''
        return self._D

    @property
    def D_D(self):
        '''N_theta x N_theta x 2 array of differences of directions:
         
          (D_D)_ijk = D_ik - D_jk
           
        Very useful in many fluxes. '''
        return self._D_D

    def dofs_on_element(self, t: int):
        '''returns the global degrees of freedom belonging to element "t"'''
        return self._T_ID_to_DOFs[t, :]
        
    def global_direction(self, n: int | int_array):
        '''computes the direction of planewave corresponding global index n'''
        m = n % self.N_theta
        return self._D[m, :]

    def plot_directions(self, ax: Axes | None = None):
        if ax is None:
            _, ax = plt.subplots()
        
        ax.plot(self._D[:, 0], self._D[:, 1], '.')
        ax.axis('square')
        plt.show()


def LinearlySpacedBasis(N_elements: int, k: float, N_theta: int, refractive_index: Elementwise[Any, float], theta_0: float = 0. ) -> PlanewaveBasis:
    """
    Construct a plane-wave basis with uniformly distributed directions.

    Parameters
    ----------
    N_elements : int
        Number of mesh elements.

    k : float
        Wavenumber.

    N_theta : int
        Number of plane-wave directions per element.

    theta_0 : float, default=0.
        Angular offset applied to all directions.

    Returns
    -------
    PlanewaveBasis
        Plane-wave basis with evenly spaced angular directions.

    Notes
    -----
    The propagation angles are generated as:

    :contentReference[oaicite:3]{index=3}

    for:

    :contentReference[oaicite:4]{index=4}
    """
    thetas = np.linspace(0, 2*np.pi, N_theta, endpoint=False) + theta_0
    return PlanewaveBasis(N_elements=N_elements, refractive_index=refractive_index, thetas=thetas, k=k)


def WaveguideTayloredBasis(N_elements: int, k: float, H: float, direction: str="right") -> PlanewaveBasis:
    """
    Construct a waveguide-adapted plane-wave basis.

    The basis directions are selected from propagating waveguide modes
    associated with a domain of height ``H``.

    Parameters
    ----------
    N_elements : int
        Number of mesh elements.

    k : float
        Wavenumber.

    H : float
        Waveguide height.

    direction : {"right", "left"}, default="right"
        Preferred propagation direction of the basis.

    Returns
    -------
    PlanewaveBasis
        Plane-wave basis adapted to propagating modal directions.

    Notes
    -----
    The modal angles are computed from the waveguide dispersion relation:

    :contentReference[oaicite:5]{index=5}

    Only propagating modes satisfying:

    :contentReference[oaicite:6]{index=6}

    are included in the basis.

    This construction is particularly useful for Trefftz discretizations
    of waveguide scattering problems.
    """

    N_P = int(k*H/np.pi)
    theta_p = np.array([np.asin(n*np.pi/(k*H)) for n in range(1, N_P+1)])  # n lamdba / H = sin(theta), theta = asin( n lambda / H) = asin( n pi / kh)
    theta_m = -theta_p
    thetas = np.concatenate([theta_m, np.array([0]), theta_p])
    if direction == "left":
        thetas = thetas + np.pi
    return PlanewaveBasis(N_elements=N_elements, thetas=thetas, k=k)



class TrefftzBasis:
    """
    Legacy Trefftz basis implementation using uniformly distributed
    plane-wave directions.

    Parameters
    ----------
    N_elements : int
        Number of mesh elements.

    N_theta : int
        Number of directions per element.

    k : float
        Wavenumber.

    Notes
    -----
    This class provides a simplified predecessor of
    :class:`PlanewaveBasis`.

    Uniform angular sampling is used for all elements.
    """
    def __init__(self, N_elements: int, N_theta: int, k: float) -> None:
        self.N_theta = N_theta
        self.k = k
        self.N_elements = N_elements
        self._N_DOF = self.N_elements * self.N_theta

        # global numbering: triangle-major order
        self.T_ID_to_DOFs = np.arange(0, self.N_DOF, dtype=np.int64).reshape(self.N_elements, self.N_theta)

        # plane-wave directions
        thetas = np.linspace(0, 2*np.pi, N_theta, endpoint=False)
        self.D = np.column_stack([np.cos(thetas), np.sin(thetas)])

    @property
    def N_DOF(self):
        return self._N_DOF

    def dofs_on_element(self, t: int):
        '''returns the global degrees of freedom belonging to element "t"'''
        return self.T_ID_to_DOFs[t, :]
        
    def global_direction(self, n: int):
        '''computes the direction of planewave corresponding global index n'''
        m = n % self.N_theta
        return self.D[m, :]
