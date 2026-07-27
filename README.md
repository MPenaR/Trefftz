# Trefftz

A Python library implementing **Discontinuous Galerkin Trefftz (DG-Trefftz)** methods for wave propagation problems. The project provides a reusable DG-Trefftz framework together with problem-specific implementations for different physical settings.

## Repository Structure

```text
trefftz/            # Core DG-Trefftz library
cases/
├── waveguide/      # Waveguide-specific implementations
└── open_space/     # Open-space scattering implementations
notebooks/          # Usage examples and tutorials
tests/              # Unit and regression tests
```

## Installation

The library requires one mesh backend. Choose **one** of the following installation options.

### Gmsh backend

```bash
pip install "trefftz[gmsh]"
```

### NGSolve backend

```bash
pip install "trefftz[ngsolve]"
```

### Both backends (optional)

```bash
pip install "trefftz[gmsh,ngsolve]"
```

### Development installation

To install the package in editable mode:

```bash
pip install -e ".[gmsh]"
```

or

```bash
pip install -e ".[ngsolve]"
```

## Current Status

The project is under active development. The core architecture is functional, but several components are still being cleaned up or migrated from an earlier implementation.

### Implemented

* Core DG-Trefftz framework.
* Waveguide kernels with exact analytical integrals.
* Extensive tests comparing the analytical waveguide kernels against numerical integration.
* Numerical implementation of the circular NtD kernels.

### Work in Progress

* Cleaning the repository to remove remnants of the previous implementation.
* Completing the analytical implementation of the circular NtD kernels. Currently, only the local contribution with the term not involving the normal is validated against the numerical implementation.
* Improving the performance of the numerical circular NtD assembly, which is currently significantly slower than the analytical waveguide implementation.
* Migrating the remaining kernels and right-hand sides from the previous codebase, including:

  * Neumann boundary condition right-hand side.
  * Dirichlet boundary condition implementation.
  * Additional kernels that are already available in the legacy implementation.

## Examples

Example notebooks demonstrating the usage of the library can be found in the `notebooks/` directory.

## Testing

Run the test suite with:

```bash
pytest
```

The tests include comparisons between analytical and numerical implementations.

