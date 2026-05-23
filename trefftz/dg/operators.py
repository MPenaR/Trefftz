from typing import Protocol
import numpy as np


class InteriorFluxOperator(Protocol):

    def assemble(self, edge, basis, k: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        ...

class BoundaryLocalFluxOperator(Protocol):

    def assemble(self, edge, basis, k: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        ...

class BoundaryGlobalFluxOperator(Protocol):

    def assemble_boundary_global(self, mesh, basis, k: float) -> sp.sparray:
        ...

