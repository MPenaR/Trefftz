# Data Structures

The mesh basically exports points, edges and triangles as numpy arrays. `points` is a simple (N,2) array of dtype np.float64, however both edges and triangles are structured arrays.


## Process

In order to allow to reusability the code is split in several classes minimizing coupling: 
- Domain: geometric definition, mesh, regions...
        - Waveguide: R, H
            - CleanWaveguide: a constructor
        - Unbounded: R
- PhysicalModel: Helmholtz, Maxwell, Lamb-waves,  sets boundary conditions, material properties, depends on a domain
- 




Sigma_L and Sigma_R should be defined separately (not as Sigma) so the assembler knows that there is no coupling among them


I need the fluxes indexed from 0 to 1 to be able to construct the numpy structured arrays BUT, once they are constructed I don't use the edge ID any more, as each edges element carries its normal its neighbours etc... so I can constructs subsets of edges as slices of edges