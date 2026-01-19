'''would be moved into dg, after making it generic but right now its waveguide-specific.'''
from trefftz.dg.basis2 import TrefftzBasis
from trefftz.numpy_types import complex_array
from .base import WaveguideRegions
from trefftz.dg.serial_fluxes2 import SoundHard, InnerEasy, Radiating_local, mode_RHS, Radiating_nonlocal
#from scipy.sparse import coo_matrix, csr_matrix, spmatrix
from trefftz.dg.fluxes import FluxType
from scipy.sparse import coo_array, csr_array, sparray
import numpy as np
from typing import Mapping

from .RHS_types import RSH_type

d_1 = 0.5
d_2 = 0.5
a = 0.5
b = 0.5



helmholtz_fluxes = {
    FluxType.SOUNDHARD: SoundHard,
    FluxType.RADIATING: Radiating_local,
    FluxType.TRANSMISSION: InnerEasy
}


flux_parameters = {
    FluxType.SOUNDHARD: {"d_1": 0.5},
    FluxType.RADIATING: {"d_2": 0.5},
    FluxType.TRANSMISSION: {"a": 0.5, "b": 0.5},
}


def SerialAssembleMatrix(edges, basis: TrefftzBasis, NtD_modes: int, boundary_conditions: Mapping[WaveguideRegions, FluxType]) -> sparray:
    N_DOF = basis.N_DOF
    I = []
    J = []
    V = []
    k = basis.k

    # LOCAL FLUXES
    # BOUNDARIES

    for edge in edges:
        match edge["region"]:
            case WaveguideRegions.GAMMA:
                flux_type = boundary_conditions[edge["region"]]
                T1, T2 = edge["triangles"]
                for m in basis.dofs_on_element(T1):
                    d_m = basis.global_direction(m)
                    for n in basis.dofs_on_element(T1):
                        d_n = basis.global_direction(n)
                        I.append(m)
                        J.append(n)
                        flux_function = helmholtz_fluxes[flux_type]
                        params = flux_parameters[flux_type]
                        V.append(flux_function(d_m=d_m, d_n=d_n, k=k, edge=edge, flux_parameters=params))

            case WaveguideRegions.SIGMA_L | WaveguideRegions.SIGMA_R:
                flux_type = boundary_conditions[edge["region"]]
                T1, T2 = edge["triangles"]
                for m in basis.dofs_on_element(T1):
                    d_m = basis.global_direction(m)
                    for n in basis.dofs_on_element(T1):
                        d_n = basis.global_direction(n)
                        I.append(m)
                        J.append(n)
                        flux_function = helmholtz_fluxes[flux_type]
                        params = flux_parameters[flux_type]
                        V.append(flux_function(d_m=d_m, d_n=d_n, k=k, edge=edge, flux_parameters=params))

# TRANSMISSION

            case _:
                flux_type = FluxType.TRANSMISSION
                T1, T2 = edge["triangles"]
                for m in basis.dofs_on_element(T1):
                    d_m = basis.global_direction(m)
                    for n in basis.dofs_on_element(T1):
                        d_n = basis.global_direction(n)
                        I.append(m)
                        J.append(n)
                        flux_function = helmholtz_fluxes[flux_type]
                        params = {"a": 0.5, "b": 0.5}
                        V.append(flux_function(d_m=d_m, d_n=d_n, k=k, edge=edge, flux_parameters=params))
                for m in basis.dofs_on_element(T1):
                    d_m = basis.global_direction(m)
                    for n in basis.dofs_on_element(T2):
                        
                        d_n = basis.global_direction(n)
                        I.append(m)
                        J.append(n)
                        flux_function = helmholtz_fluxes[flux_type]
                        params = {"a": -0.5, "b": -0.5}
                        V.append(flux_function(d_m=d_m, d_n=d_n, k=k, edge=edge, flux_parameters=params))
                for m in basis.dofs_on_element(T2):
                    d_m = basis.global_direction(m)
                    for n in basis.dofs_on_element(T1):
                        d_n = basis.global_direction(n)
                        I.append(m)
                        J.append(n)
                        flux_function = helmholtz_fluxes[flux_type]
                        params = {"a": 0.5, "b": 0.5}
                        V.append(-flux_function(d_m=d_m, d_n=d_n, k=k, edge=edge, flux_parameters=params))
                for m in basis.dofs_on_element(T2):
                    d_m = basis.global_direction(m)
                    for n in basis.dofs_on_element(T2):
                        d_n = basis.global_direction(n)
                        I.append(m)
                        J.append(n)
                        flux_function = helmholtz_fluxes[flux_type]
                        params = {"a": -0.5, "b": -0.5}
                        V.append(-flux_function(d_m=d_m, d_n=d_n, k=k, edge=edge, flux_parameters=params))


# NON LOCAL
    for edge_u in edges[np.where(edges["region"]==WaveguideRegions.SIGMA_L)[0]]:
            flux_type = boundary_conditions[edge_u["region"]]
            T_u, _ = edge_u["triangles"]
            for n in basis.dofs_on_element(T_u):
                d_n = basis.global_direction(n)
                for edge_v in edges[np.where(edges["region"]==WaveguideRegions.SIGMA_L)[0]]:
                    T_v, _ = edge_v["triangles"]
                    for m in basis.dofs_on_element(T_v):
                        d_m = basis.global_direction(m)
                        I.append(m)
                        J.append(n)
                        params = flux_parameters[flux_type]
                        V.append(Radiating_nonlocal(d_m=d_m, d_n=d_n, edge_v=edge_v, edge_u=edge_u, k=k, NtD_modes=NtD_modes, flux_parameters=params))


    A = csr_array(coo_array((V, (I, J)), shape=(N_DOF, N_DOF)))

    return A

def SerialAssembleRHS(edges, basis: TrefftzBasis, NtD_modes: int, RHS: RSH_type, RHS_params: Mapping[str, int | float]) ->sparray:
    N_DOF = basis.N_DOF
    I = []
    J = []
    V = []
    k = basis.k

    match RHS:
        case RSH_type.MODE:
            t = RHS_params["t"]
            I = []
            V = []
            for edge in edges[np.where(edges["region"]==WaveguideRegions.SIGMA_L)[0]]:
                T, _ = edge["triangles"]
                for m in basis.dofs_on_element(T):
                    d_m = basis.global_direction(m)
                    I.append(m)
                    V.append(mode_RHS(d_m=d_m, edge=edge, k=k, H=1., d_2=0.5, t=t))

    b = coo_array((V, (I,)), shape=(N_DOF,))

    return b



# def AssembleMatrix(V : TrefftzSpace,  Edges : tuple[Edge], 
#                    H : float, a = 0.5, b = 0.5, d_1 = 0.5, d_2 = 0.5, 
#                    Np=10, full_matrix = False) -> spmatrix:
#     '''Assembles de matrix for the bilinear form.
#     a, b, d_1 and d_2 are the coefficients of the regularizing terms.
#     Use full_matrix = Truee if the returned matrix should NOT be a sparse
#     matrix.'''


#     N_DOF = V.N_DOF

#     values = []
#     i_index = []
#     j_index = []


#     Side_edges  : dict[EdgeType, list] = { EdgeType.SIGMA_L : [], EdgeType.SIGMA_R : []} 
    
#     for E in Edges:
#         match E.Type:
#             case EdgeType.SIGMA_L | EdgeType.SIGMA_R:
#                 Side_edges[E.Type].append(E)
#             case _:
#                 pass

#     N_Edges = len(Edges)

    
#     # filling the vectors
#     if np.isscalar(a):
#         a_vec = np.full(N_Edges,a)
#     else:
#         a_vec = a 
        
#     if np.isscalar(b):
#         b_vec = np.full(N_Edges,b)
#     else:
#         b_vec = b 
    
#     if np.isscalar(d_1):
#         d_1_vec = np.full(N_Edges,d_1)
#     else:
#         d_1_vec = d_1 
    

#     if np.isscalar(d_2):
#         d_2_vec = np.full(N_Edges,d_2)
#     else:
#         d_2_vec = d_2
    


#     Phi = V.TrialFunctions
#     Psi = V.TestFunctions # currently the same spaces 

#     k = V.kappa

#     for (E, a, b, d_1, d_2) in zip(Edges, a_vec, b_vec, d_1_vec, d_2_vec):
#         match E.Type:
#             case EdgeType.INNER:
#                 K_plus, K_minus = E.Triangles

#                 for n in V.DOF_range[K_plus]:
#                     phi = Phi[n]
#                     for m in V.DOF_range[K_plus]:
#                         psi = Psi[m]
#                         i_index.append(m)
#                         j_index.append(n)
#                         values.append(Inner_term_general(phi, psi, E, k, a, b))

#                 for n in V.DOF_range[K_minus]:
#                     phi = Phi[n]
#                     for m in V.DOF_range[K_plus]:
#                         psi = Psi[m]
#                         i_index.append(m)
#                         j_index.append(n)
#                         values.append(Inner_term_general(phi, psi, E, k, -a, -b))


#                 for n in V.DOF_range[K_plus]:
#                     phi = Phi[n]
#                     for m in V.DOF_range[K_minus]:
#                         psi = Psi[m]
#                         i_index.append(m)
#                         j_index.append(n)
#                         values.append(-Inner_term_general(phi, psi, E, k, a, b))

#                 for n in V.DOF_range[K_minus]:
#                     phi = Phi[n]
#                     for m in V.DOF_range[K_minus]:
#                         psi = Psi[m]
#                         i_index.append(m)
#                         j_index.append(n)
#                         values.append(-Inner_term_general(phi, psi, E, k, -a, -b))


#             case EdgeType.GAMMA:
#                 K = E.Triangles[0]
#                 for m in V.DOF_range[K]:
#                     psi = Psi[m]
#                     for n in V.DOF_range[K]:
#                         phi = Phi[n]
#                         i_index.append(m)
#                         j_index.append(n)
#                         values.append(Gamma_term(phi, psi, k, E, d_1))
                    

#             case EdgeType.D_OMEGA | EdgeType.COVER:
#                 K = E.Triangles[0]
#                 for m in V.DOF_range[K]:
#                     psi = Psi[m]
#                     for n in V.DOF_range[K]:
#                         phi = Phi[n]
#                         i_index.append(m)
#                         j_index.append(n)
#                         values.append(sound_soft_term(phi, psi, k, E, a))


#             case EdgeType.SIGMA_L | EdgeType.SIGMA_R:
#                 K = E.Triangles[0]
#                 for n in V.DOF_range[K]:
#                     phi = Phi[n]
#                     for E_other in Side_edges[E.Type]:
#                         K_other = E_other.Triangles[0]
#                         if E_other == E:
#                             for m in V.DOF_range[K_other]:
#                                 psi = Psi[m]
#                                 i_index.append(m)
#                                 j_index.append(n)
#                                 S = Sigma_local(phi, psi, k, E, d_2) + Sigma_nonlocal(phi, psi, E, E, k, H, d_2, Np=Np)
#                                 values.append(S)
#                         else:
#                             for m in V.DOF_range[K_other]:
#                                 psi = Psi[m]
#                                 i_index.append(m)
#                                 j_index.append(n)
#                                 values.append(Sigma_nonlocal(phi, psi, E, E_other, k, H, d_2, Np=Np))
#     if V.absorbing:                    
#         for T in V.ScattererTriangles:
#             r_A, r_B, r_C = T.A, T.B, T.C
#             T_index = T.index
#             for n in V.DOF_range[T_index]:
#                 phi = Phi[n]
#                 for m in V.DOF_range[T_index]:
#                     psi = Psi[m]
#                     i_index.append(m)
#                     j_index.append(n)
#                     values.append(absorption_term( phi, psi, r_A, r_B, r_C, k))

          
#     A = coo_matrix( (values, (i_index, j_index)), shape=(N_DOF,N_DOF))
#     A = csr_matrix(A)

#     if full_matrix:
#         A = A.toarray()

#     return A

