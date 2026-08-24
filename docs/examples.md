# Running examples

Run examples from the repository root with the `toast` Conda environment active and the package installed in editable mode.

## Layered wall

```bash
python examples/layered_wall.py
```

This is the closest example to the original prototype. It demonstrates:

- `LayeredWallSpec` and the conforming layered Gmsh builder;
- semantic layer-to-material mapping;
- exterior/interior convection boundaries;
- transient heat conduction;
- VTX/ADIOS2 output;
- material, temperature, and profile visualization.

Use this example when checking changes to the transient solver or layered-wall mesh generation.

## Manufactured 2-D solution

```bash
python examples/manufactured.py
```

This example solves a homogeneous 2-D problem with a source selected so that an exact sinusoidal temperature field is known. It is the visualization counterpart to the manufactured-solution convergence test.

The manufactured solution is useful because it exercises genuine two-dimensional gradients and curvature without involving material interfaces or complex geometry.

## Wall/floor/foundation junction

```bash
python examples/wall_floor_foundation.py
```

This is the current complex geometry example. It includes:

- continuous exterior insulation from the wall down the foundation;
- sheathing, cavity insulation, service chase, and drywall above the floor;
- a wood structural floor/rim region extending into the interior;
- a concrete foundation beneath the floor;
- exterior-temperature conditions on the outer insulation face and exposed below-floor surfaces;
- interior-temperature conditions on the drywall face and floor top;
- unclassified artificial truncation boundaries that remain adiabatic.

The geometry uses `RectangleRegion` and `BoundaryRule`, then calls the same `solve_steady_heat()` used by the simpler benchmark cases.

The example prints:

- global cell count;
- wall/foundation dimensions;
- exterior and interior integrated heat flows;
- relative energy imbalance;
- 2-D thermal coupling coefficient.

It also writes `wall_floor_foundation.msh` and, in serial runs, displays material and temperature fields.

## Stud-wall benchmark

There is no standalone example script yet, but `tests/test_stud_wall_benchmark.py` is a useful executable reference for a two-material thermal bridge. It compares the FEM effective U-value against all-insulation, all-wood, and area-weighted parallel-path reference values.

## Running with MPI

The solver and mesh conversion APIs carry MPI communicators, and the examples use `MPI.COMM_WORLD`. Numerical work can be launched under MPI where supported, for example:

```bash
mpirun -n 2 python examples/wall_floor_foundation.py
```

Interactive visualization is currently guarded to serial runs in the junction example. Point-sampling helpers are also primarily intended for serial post-processing.
