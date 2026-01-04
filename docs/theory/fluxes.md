# Fluxes

All along we will be using $l$ as the length of the edge going from $\mathbf{P}$ to $\mathbf{Q}$. Its mid point will be denoted as $\mathbf{M}$ and the tangent and normal as $\boldsymbol{\tau}$ and $\mathbf{n}$ respectively.

The trial and test functions will be plane waves of the form:

```{math}
e^{ik\mathbf{d}\cdot\mathbf{x}}
```

## Sound-hard boundary condition

If the edge belongs to a boundary where we have imposed a sound-hard boundary condition, i.e.: 

```{math}
\nabla u \cdot \mathbf{n} = 0
```

then the corresponding flux becomes:

```{math}
\int_E \left(u + \frac{\mathfrak{d}_1}{ik} \nabla u \cdot \mathbf{n}_E\right)\overline{\nabla v \cdot \mathbf{n}_E} \,\mathrm d \ell
```

So the term corresponding to test function $\psi_m$ and trial function $\varphi_n$ is:

```{math}
\int_E \left(\varphi_n(\mathbf{x}) + \frac{\mathfrak{d}_1}{ik} \nabla \varphi_n(\mathbf{x}) \cdot \mathbf{n}_E\right)\overline{\nabla \psi_m(\mathbf{x}) \cdot \mathbf{n}_E} \,\mathrm d \ell = -ik\mathbf{d}_m \cdot \mathbf{n}_E \left(1 + \mathfrak{d}_1\mathbf{d}_n \cdot \mathbf{n}_E  )\int_E \varphi_n(\mathbf{x} \right)\overline{ \psi_m(\mathbf{x}) } \,\mathrm d \ell 
```
