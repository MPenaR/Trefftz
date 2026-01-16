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
