# Introduction

`trefftz` is a Python library implementing a discontinuous Galerkin Trefftz method for the numerical solution of the Helmholtz equation.

The project separates the general numerical machinery from problem-specific components. The `trefftz` package contains the reusable components required to discretize and solve Helmholtz problems, while the `cases` directory contains implementations specific to particular physical configurations, currently open-space and waveguide problems.

This separation allows the same numerical machinery to be reused across different physical configurations and makes it possible to test general components independently of any particular application.