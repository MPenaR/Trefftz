'''module for defining the problem  base class'''
from typing import TYPE_CHECKING, Optional, Callable, TypeVar, Generic

if TYPE_CHECKING:
    from trefftz.dg.functions import ComplexFunction  # consider changing it to "field"

from trefftz.mesh import TrefftzMesh


from trefftz.numpy_types import float_array, complex_array
from trefftz.dg.fluxes import FluxType
from trefftz.dg.basis import TrefftzBasis

from enum import IntEnum
from abc import abstractmethod
import numpy as np
from typing import Protocol

R = TypeVar("R", bound=IntEnum)


class Domain(Protocol[R]):
    mesh: TrefftzMesh
    regions: type[R]

    def __init__(self, mesh: TrefftzMesh):
        self.mesh = mesh
    
    # @abstractmethod
    # def generate_Trefftzmesh(self):
    #     """Produce TrefftzMesh."""


D = TypeVar("D")


class Problem:
    '''Class which manages defining domains and setting boundary conditions'''

    def __init__(self, mesh: "TrefftzMesh", boundary_conditions_map: Optional[dict[IntEnum, FluxType]] | None,
                  k: float, verbose: bool = False):
        self.domain = mesh
        self.verbose = verbose
        if boundary_conditions_map:
            self.set_boundary_conditions(boundary_conditions_map)
        self.k = k

    @abstractmethod
    def build_geometry(self):
        """Create gmsh geometry or analytic geometry."""

    @abstractmethod
    def generate_Trefftzmesh(self):
        """Produce TrefftzMesh."""

    def set_boundary_conditions(self, boundary_conditions_map: dict[IntEnum, FluxType]):
        self.boundary_conditions_map = boundary_conditions_map
        """Assign BCs into BoundaryConditionModel."""
        for region, flux in boundary_conditions_map.items():
            self.domain.edges["flux_type"][self.domain.edges["region"] == region] = flux
        if np.all(self.domain.edges["flux_type"] >= 0):
            if self.verbose: 
                print('Problem ready for assembly')
            self.domain.ready_for_assemble = True

    @abstractmethod
    def plot(self, figsize: Optional[tuple[int, int]] = (16, 2), line_width: Optional[int] = 1):
        """Visualize geometry/mesh/regions."""

    @abstractmethod
    def plot_field(self, u: Callable[[float_array, float_array], complex_array],
                   N: int = 100, figsize: Optional[tuple[int, int]] = (16, 2), real_part: bool = False):
        """Visualize fields defined on the mesh."""

    def assemble(self, basis):
        pass


class PhysicalModel:
    def __init__(self, domain: Domain, boundary_conditions, materials):
        self.domain = domain
        self.boundary_conditions = boundary_conditions
        self.materials = materials
    pass


class AbstractProblem(Generic[D, R]):
    '''Class which manages defining domains and setting boundary conditions'''
    domain: D
    physics: PhysicalModel
    verbose: bool
    regions: R
    basis: TrefftzBasis
    k: float


    def __init__(self, domain: D, model: PhysicalModel, basis: TrefftzBasis, regions: R, k: float, verbose: bool = False):
        self.domain = domain
        self.physics = model  # model(domain,k, materials: dict)
        self.verbose = verbose
        self.regions = regions
        self.basis = basis
        self.k = k

    @abstractmethod
    def plot(self, figsize: Optional[tuple[int, int]] = (16, 2), line_width: Optional[int] = 1):
        """Visualize geometry/mesh/regions."""

    @abstractmethod
    def plot_field(self, u: Callable[[float_array, float_array], complex_array],
                   N: int = 100, figsize: Optional[tuple[int, int]] = (16, 2), real_part: bool = False):
        """Visualize fields defined on the mesh."""

    def assemble(self):
        pass

    # def set_boundary_conditions(self, boundary_conditions_map: dict[IntEnum, FluxType]):
    #     self.boundary_conditions_map = boundary_conditions_map
    #     """Assign BCs into BoundaryConditionModel."""
    #     for region, flux in boundary_conditions_map.items():
    #         self.domain.edges["flux_type"][self.domain.edges["region"] == region] = flux
    #     if np.all(self.domain.edges["flux_type"] >= 0):
    #         if self.verbose:
    #             print('Problem ready for assembly')
    #         self.domain.ready_for_assemble = True
