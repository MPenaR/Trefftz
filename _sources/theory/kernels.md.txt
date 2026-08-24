# Kernels

This page describes the mathematical expressions used to evaluate the edge
integrals appearing in the Trefftz-DG formulation.

The kernels considered in the library are

- `I_uv`: the integral of the product of two basis functions.
- `I_duv`: the integral involving the normal derivative of one of the
  basis functions.

Both kernels are defined for the two types of edges supported by the library:

- straight **segments**;
- **circular arcs**.

The resulting expressions are summarized below.

## Kernel notation

Introduce here the general notation used throughout this page, including:

- the two functions $u=e^{ik\mathbf{d}_u\cdot\mathbf{x}}$ and $v=e^{ik\mathbf{d}_v\cdot\mathbf{x}}$;
- the edge $E$ where the functions are supported;
- the  normal $\mathbf{n}$;
- the midpoint of the edge $\mathbf{M}$
- the tangent vector $\boldsymbol{\tau}$.
- the length $l$.

The kernels are defined as

### `I_uv`

$$
I_{uv} := \int_E u(\mathbf{x})\overline{v(\mathbf{x})}\,\mathrm{d}\ell
$$

### `I_duv`

$$
I_{duv} := \int_E \nabla u(\mathbf{x})\cdot{\mathbf{n}}\overline{v(\mathbf{x})}\,\mathrm{d}\ell
$$

---

## Kernels on straight segments

For a straight segment, the integrals can be evaluated analytically and both of them are indeed proportional.

### `I_uv`

$$
I_{uv} = le^{ik\left(\mathbf{d}_u-\mathbf{d}_v\right)\cdot\mathbf{M}}\mathrm{sinc}\left(\frac{kl}{2\pi}\left(\mathbf{d}_u-\mathbf{d}_v\right)\cdot\boldsymbol{\tau}\right)
$$

where $\mathrm{sinc}$ stands for the _sinus cardinalis_ function: 

$$
\mathrm{sinc}(x) = \frac{\sin(\pi x)}{\pi x}
$$

### `I_duv`

$$
I_{duv} = ik\mathbf{d}_u\cdot\mathbf{n}I_{uv}
$$

Just written for analogy with the arcs.

---

## Kernels on circular arcs

For a circular arc, the corresponding integrals are evaluated using the
series representation of a plane wave given by the Jacobi-Angers expansion.

### `I_uv`

$$
I_{uv}=R\left(J_{0}\left(kRD_{uv}\right)\left(\theta_{2}-\theta_{1}\right)+2\sum_{t=1}^{\infty}\frac{i^{t}}{t}J_{t}\left(kRD_{uv}\right)\left(\sin\left(t\left(\theta_{2}-\phi_{uv}\right)\right)-\sin\left(t\left(\theta_{1}-\phi_{uv}\right)\right)\right)\right)
$$

where $J_v$ stands for the $v$ order Bessel function of the first kind, $R$ is the radius of the arc and $D_{uv}$ and $\phi_{uv}$ are defined such that:

$$
\mathbf{d}_u - \mathbf{d}_v = D_{uv}\left(\cos(\phi_{uv})\mathbf{i} + \sin(\phi_{uv})\mathbf{j}\right)
$$

### `I_duv`

$$
I_{duv} = I(\theta_2) - I(\theta_1)
$$

where the primitive $I$ is

$$
I(\theta)=kR\left(-J_{1}\left(kRD_{uv}\right)\cos\left(\phi_{uv}-\phi_{u}\right)\theta+\sum_{p=1}^{\infty}\frac{i^{p}}{p}\left(J_{p-1}\left(kRD_{uv}\right)\sin\left(p\left(\theta-\phi_{uv}\right)+\left(\phi_{uv}-\phi_{u}\right)\right)-J_{p+1}\left(kRD_{uv}\right)\sin\left(p\left(\theta-\phi_{uv}\right)-\left(\phi_{uv}-\phi_{u}\right)\right)\right)\right)
$$

and $\phi_u$ is defined such that $\mathbf{d}_u = \left(\cos(\phi_{u})\mathbf{i} + \sin(\phi_{u})\mathbf{j}\right)$.

---

## Summary

The four kernel/geometry combinations implemented in the library are:

| Geometry | Kernel | Evaluation |
|---|---|---|
| Segment | $I_{uv}$ | Exact |
| Segment | $I_{duv}$ | Exact |
| Circular arc | $I_{uv}$ | Truncated series |
| Circular arc | $I_{duv}$ | Truncated series |

The mathematical expressions described here form the basis for the generic
kernel implementations in `trefftz`. Their numerical implementations are
tested independently against numerical integration; see
[Testing](../design/testing.md) for details.