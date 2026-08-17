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

For a straight segment, the integrals can be evaluated analytically.

### `I_uv`

$$
I_{uv} = le^{ik\left(\mathbf{d}_u-\mathbf{d}_v\right)\cdot\mathbf{M}}\mathrm{sinc}\left(\frac{kl}{2\pi}\left(\mathbf{d}_u-\mathbf{d}_v\right)\cdot\boldsymbol{\tau}\right)
$$

[Brief explanation of the notation appearing in the formula.]

### `I_duv`

$$
I_{duv} = ik\mathbf{d}_u\cdot\mathbf{n}I_{uv}
$$

[Brief explanation of the notation appearing in the formula.]

---

## Kernels on circular arcs

For a circular arc, the corresponding integrals are evaluated using the
series representation derived for the circular geometry.

### $I_{uv}$

[Latex formula for $I_{uv}$ on a circular arc]

[Brief explanation of the notation and the series.]

### $I_{duv}$

[Latex formula for $I_{duv}$ on a circular arc]

[Brief explanation of the notation and the series.]

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