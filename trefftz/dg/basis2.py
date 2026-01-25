import numpy as np
from trefftz.numpy_types import float_array, int_array
from dataclasses import dataclass
from matplotlib.axes import Axes
import matplotlib.pyplot as plt

@dataclass(slots=True)
class PlanewaveBasis:
    '''core basis class'''

    k: float
    _N_theta: int
    _N_DOF: int
    _N_elements: int
    _T_ID_to_DOFs: int_array
    _D: float_array
        
    def __init__(self, N_elements: int, thetas: float_array, k: float) -> None:
        self._N_theta = len(thetas)
        self.k = k
        self._N_elements = N_elements
        self._N_DOF = self.N_elements * self.N_theta

        # global numbering: triangle-major order
        self._T_ID_to_DOFs = np.arange(0, self.N_DOF, dtype=np.int64).reshape(self.N_elements, self.N_theta)

        # plane-wave directions
        self._D = np.column_stack([np.cos(thetas), np.sin(thetas)])

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
    

    def dofs_on_element(self, t: int):
        '''returns the global degrees of freedom belonging to element "t"'''
        return self._T_ID_to_DOFs[t, :]
        
    def global_direction(self, n: int):
        '''computes the direction of planewave corresponding global index n'''
        m = n % self.N_theta
        return self._D[m, :]

    def plot_directions(self, ax: Axes | None = None):
        if ax is None:
            _, ax = plt.subplots()
        
        ax.plot(self._D[:, 0], self._D[:, 1], '.')
        ax.axis('square')
        plt.show()


def LinearlySpacedBasis(N_elements: int, k: float, N_theta: int, theta_0: float = 0. ) -> PlanewaveBasis:
    thetas = np.linspace(0, 2*np.pi, N_theta, endpoint=False) + theta_0
    return PlanewaveBasis(N_elements=N_elements, thetas=thetas, k=k)


def TayloredBasis(N_elements: int, k: float, H: float, direction: str="right") -> PlanewaveBasis:

    N_P = int(k*H/np.pi)
    theta_p = np.array([np.asin(n*np.pi/(k*H)) for n in range(1, N_P+1)])  # n lamdba / H = sin(theta), theta = asin( n lambda / H) = asin( n pi / kh)
    theta_m = -theta_p
    thetas = np.concatenate([theta_m, np.array([0]), theta_p])
    if direction == "left":
        thetas = thetas + np.pi
    return PlanewaveBasis(N_elements=N_elements, thetas=thetas, k=k)



class TrefftzBasis:
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
    


# class CheatyBasis:
#     def __init__(self, N_elements: int, H: float, k: float) -> None:
#         '''creates a basis taylored to contain the propagatin modes of a waveguide'''
        
#         N_P = int(k*H/np.pi)
#         theta_p = np.array([np.asin(n*np.pi/(k*H)) for n in range(1,N_P+1)])     #  n lamdba / H = sin(theta), theta = asin( n lambda / H) = asin( n pi / kh)
#         theta_m = -theta_p
#         thetas = np.concatenate([theta_m, np.array([0]), theta_p])
#         self.k = k
#         self.N_elements = N_elements

#         self.D = np.column_stack([np.cos(thetas), np.sin(thetas)])
#         # self.D = np.concatenate([self.D, -self.D], axis=0)
#         self.N_theta = self.D.shape[0]
#         self._N_DOF = self.N_elements * self.N_theta



#         # global numbering: triangle-major order
#         self.T_ID_to_DOFs = np.arange(0, self.N_DOF, dtype=np.int64).reshape(self.N_elements, self.N_theta)

#         # plane-wave directions
 
   
#     @property
#     def N_DOF(self):
#         return self._N_DOF

#     def dofs_on_element(self, t: int):
#         return self.T_ID_to_DOFs[t, :]
        
#     def global_direction(self, n: int):
#         '''computes the direction of planewave corresponding global index n'''
#         m = n % self.N_theta
#         return self.D[m, :]

