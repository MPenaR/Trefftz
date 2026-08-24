# Problem object

The problem class is an abstract class or template (we will see), so I can hace more specific problems like:

- Helmholtz waveguide problem
- Helmholtz free-space problem
- Lamb plate problem

All of them share some features:

- A domain needs to be defined
- A function basis needs to be chosen
- Some physical model, with is boundary conditions are set
- A matrix A and a vector b are assembled
- The linear system Ax=b is solved
- The solution is compared to an exact one,
- Error metrics are computed

As such this class should be a top class that calls specific domain creators, basis, assemblers etc...

This is interesting because for example, the exact solutions are something which is specific to the pair (physics + domain). The fundamental solution for the helmholtz problem in a waveguide is not the same than in free space, nor is the same for the Helmholtz equation or the Maxwell equations.

## Waveguide problem

A waveguide is a class in itself. Why? because it has parameters defining it ($R$ and $H$), it has specific plotting functions and it has modes (exact solutions). We do not have an analogous class for the free-space problem, we dont have a "free-space". Perhaps we should, if we say that the truncation distance $R$ defines the waveguide, then $R$ would also define the truncation of the free-space (i.e. a centered circle).

In any case, I don't think "problem" should expect to know if its working with a waveguide or in free-space. Problem knows about:

- The mesh.
- The basis (or its directions)
- Physics: the boundary conditions imposed (i.e. a dictionary between known regions and boundary condition type)
- Numerics: the mapping between boundary condition and kernel ( specially important to set the correct NtD kernels, circular/waveguide)
- Assembling A and b
- Solving the system 
- Computing errors against a given solution (perhaps an exact one). 
