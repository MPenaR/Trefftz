'''Defining an managing the Trefftz basis (by now only plane waves, fixed P)'''

import numpy as np
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from trefftz.mesh import TrefftzMesh


class TrefftzBasis:
    def __init__(self, mesh: "TrefftzMesh", N_theta: int, k: float) -> None:
        self.mesh = mesh
        self.N_theta = N_theta
        self.k = k

        self.N_triangles = mesh.n_triangles
        self._N_DOF = self.N_triangles * self.N_theta

        # global numbering: triangle-major order
        self.T_ID_to_DOFs = np.arange(0, self.N_DOF, dtype=np.int64).reshape(self.N_triangles, self.N_theta)

        # plane-wave directions
        thetas = np.linspace(0, 2*np.pi, N_theta, endpoint=False)
        self.D = np.column_stack([np.cos(thetas), np.sin(thetas)])

    @property
    def N_DOF(self):
        return self._N_DOF

    def dofs_on_triangle(self, t: int):
        return self.T_ID_to_DOFs[t, :]

    # @property
    # def directions(self):
    #     return self.D
