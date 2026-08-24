# Verification and testing

The test suite is structured so that increasingly complex parts of the software stack are introduced one at a time.

Run everything with:

```bash
pytest -q
```

For printed convergence/benchmark diagnostics:

```bash
pytest -v -s
```

## Analytic layered-wall test

`tests/test_analytic.py` verifies the pure-Python layered resistance, heat-flux, and temperature-profile calculations independently of DOLFINx.

This catches sign, resistance, and material-mapping mistakes without involving the FEM stack.

## Uniform-equilibrium transient invariant

`tests/test_uniform_equilibrium.py` initializes a single-material wall and both convection environments to the same temperature. The exact solution is constant for all time.

The test checks that the backward-Euler transient implementation preserves that equilibrium to a tight numerical tolerance.

## Steady layered wall versus analytic solution

`tests/test_steady_layered_wall.py` solves a multi-layer wall with convection boundaries and compares a mid-height numerical profile with the exact one-dimensional piecewise-linear solution.

This simultaneously checks:

- tagged material-property fields;
- convection weak-form terms;
- steady solve;
- point sampling;
- agreement with the series-resistance model.

## 2-D manufactured solution

`tests/test_manufactured_2d.py` uses an exact solution of the form

```text
T(x,y) = T0 + A sin(pi x/Lx) sin(pi y/Ly)
```

with a matching volumetric source. For P1 elements the L2 error should approach second-order convergence, so halving the mesh spacing should reduce the error by approximately a factor of four.

The current test compares coarse and fine meshes and requires an error ratio greater than 3.

## Patchwork consistency test

`tests/test_patchwork_consistency.py` subdivides a rectangle into four independently tagged Gmsh regions but assigns the same conductivity to all four. It repeats the manufactured-solution convergence check.

This isolates the patchwork mesh/tag path:

```text
Gmsh rectangles -> fragmentation -> cell tags -> DG0 k -> FEM solve
```

without introducing a real material discontinuity.

## Stud-wall thermal-bridge benchmark

`tests/test_stud_wall_benchmark.py` introduces a wood stud and insulation with convection boundaries. It checks:

- steady global energy conservation using Robin boundary heat flow;
- FEM U-value greater than an all-insulation wall;
- FEM U-value less than an all-wood wall;
- reasonable agreement with a parallel-path engineering approximation.

This is a physical regression benchmark rather than an exact analytic solution.

## What should be tested next

For more complex junction geometries, add two convergence studies before relying on reported thermal-bridge metrics:

1. **mesh convergence**: reduce local/global element size and confirm heat-flow metrics stabilize;
2. **domain-size convergence**: increase artificial wall height, foundation depth, and floor extension until junction heat flow is insensitive to truncation boundaries.

When new geometry functionality is added, prefer a test that isolates geometry/tag behavior before coupling it to new physics.
