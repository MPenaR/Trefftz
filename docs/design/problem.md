# Problem object

The problem cass is an abstract class or template (we will see), so I can hace more specific problems like:

- Helmholtz waveguide problem
- Helmholtz unbounded problem
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