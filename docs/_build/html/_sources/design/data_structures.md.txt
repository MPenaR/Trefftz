# Data Structures

The mesh basically exports points, edges and triangles as numpy arrays. `points` is a simple (N,2) array of dtype np.float64, however both edges and triangles are structured arrays.

## Mesh

- A point $\mathbf{x} \in \mathbb{R}^n$ is _never_ in a subset with 0 $n$-Lebesgue measure. Hence `get_ID` only returns triangles.

- A mesh has two types of structures/philosophies:
  - For mesh creation and manipulation we use pythonic datacontainers and methods ( collections, lists comprehensions, dictionaries...). They are not vectorized but they are only called a couple of times during initialization of the mesh.
  - Mesh only exports numpy structured arrays for fast manipulation by the DG-Trefftz assembler. These are vectorized, the assembler knows nothing about the mesh internals.

<!-- Mesh internals: 
- an edge IS a couple of nodes, they could be a frozenset, they could be hashable.

A mesh should export two different sets of edges: inner edges and outer edges
The former belong to two triangles while the last only to one. Furthermore the this two subset of edges should be further separated by boundary conditions (by boundary definition) -->

## Process

In order to allow to reusability the code is split three packages:

- trefftz: general functions, independent of the particular problem.
- cases:
  - waveguide: functions only meaningfull for a waveguide problem: waveguide modes, waveguide plots, NtD...
  - open_space: functions only meaningfull for an open space problem: Fourier modes, circular plots, NtD...

<!-- - Domain: geometric definition, mesh, regions...
        - Waveguide: R, H
            - CleanWaveguide: a constructor
        - Unbounded: R
- PhysicalModel: Helmholtz, Maxwell, Lamb-waves,  sets boundary conditions, material properties, depends on a domain -->
<!-- - 

Sigma_L and Sigma_R should be defined separately (not as Sigma) so the assembler knows that there is no coupling among them


I need the fluxes indexed from 0 to 1 to be able to construct the numpy structured arrays BUT, once they are constructed I don't use the edge ID any more, as each edges element carries its normal its neighbours etc... so I can constructs subsets of edges as slices of edges -->