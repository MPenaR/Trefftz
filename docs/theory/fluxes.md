# Fluxes

All along we will be using $l_E$ as the length of the edge going from $\mathbf{P}$ to $\mathbf{Q}$. Its mid point will be denoted as $\mathbf{M}_E$ and the tangent and normal as $\boldsymbol{\tau}_E$ and $\mathbf{n}_E$ respectively.

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
\int_E \left(u + \frac{\mathfrak{d}_1}{ik} \nabla u \cdot \mathbf{n}_E\right)\overline{\nabla v \cdot \mathbf{n}_E} \,\mathrm{d} \ell
```

So the term corresponding to test function $\psi_m$ and trial function $\varphi_n$ is:

```{math}
\begin{align*}
\int_E \left(\varphi_n(\mathbf{x}) + \frac{\mathfrak{d}_1}{ik} \nabla \varphi_n(\mathbf{x}) \cdot \mathbf{n}_E\right)\overline{\nabla \psi_m(\mathbf{x}) \cdot \mathbf{n}_E} \,\mathrm d \ell &= -ik\mathbf{d}_m \cdot \mathbf{n}_E \left(1 + \mathfrak{d}_1\mathbf{d}_n \cdot \mathbf{n}_E \right)\int_E \varphi_n(\mathbf{x})\overline{ \psi_m(\mathbf{x}) } \,\mathrm {d} \ell \\
&= -ik\mathbf{d}_m \cdot \mathbf{n}_E \left(1 + \mathfrak{d}_1\mathbf{d}_n \cdot \mathbf{n}_E \right)\int_E e^{ik\left(\mathbf{d}_n-\mathbf{d}_m\right)\cdot\mathbf{x}} \,\mathrm d \ell \\
&= -ikl_E\mathbf{d}_m \cdot \mathbf{n}_E \left(1 + \mathfrak{d}_1\mathbf{d}_n \cdot \mathbf{n}_E \right) e^{ik\left(\mathbf{d}_n-\mathbf{d}_m\right)\cdot\mathbf{M}}  \mathrm{sinc}\left(\frac{kl}{2\pi}\left(\mathbf{d}_n-\mathbf{d}_m\right)\cdot\boldsymbol{\tau}_E\right) 
\end{align*}
```

This clearly shows that each edge with a sound-hard boundary condition contributes to the whole assembled matrix with a full $N_\theta \time N_\theta$ block in the main diagonal (corresponding to the triangle to the edge belongs) whose shape can be computed in a vectorized manner as:

```{math}
G = ikl_E\mathbf{d}_m \cdot \mathbf{n}_E \left(1 + \mathfrak{d}_1\mathbf{d}_n \cdot \mathbf{n}_E \right) e^{-ik\Delta D\mathbf{M}}\circ \mathrm{sinc}\left(\frac{kl}{2\pi}(\Delta D\mathbf{\tau}_E)\right)
```

where $\Delta D$ is the  $N_\theta \times N_\theta \times 2 $ array constructed as the "outer difference" of the set of directions:

```{math}
\mathscr{D}_{mnl} := D_{ml} - D_{nl}
```

$D$ beint the $N_\theta\times 2$ array of directions:

```{math}
D = \begin{bmatrix}
\cos(\theta_0) & \sin(\theta_0) \\
\cos(\theta_1) & \sin(\theta_1) \\
\vdots & \vdots \\
\cos(\theta_{N_\theta-1}) & \sin(\theta_{N_\theta-1}) \\
\end{bmatrix}
```

The product with $\mathbf{n}$ or $\mathbf{M}$ is computed along the last axis and $\circ$ denotes the Haddamard or element-wise product of two same matrices with the same shape.