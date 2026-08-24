# Thermal Object Analysis and Simulation Toolkit (TOAST)

`TOAST` is a small Python package for two-dimensional finite-element heat-transfer models of building-envelope assemblies and thermal bridges. It combines **Gmsh** geometry/meshing with **DOLFINx/FEniCSx** heat-equation solvers, while keeping geometry, material properties, boundary conditions, numerical solution, and post-processing as separate concerns.

The current package supports:

- steady and transient heterogeneous heat conduction in 2-D;
- convection (Robin) boundary conditions;
- layered-wall meshes with conforming material interfaces;
- rectangular patchwork geometries with Gmsh fragmentation;
- geometric boundary classification with reusable `BoundaryRule` predicates;
- material-property assignment through DOLFINx cell tags;
- integrated heat-flow, energy-balance, U-value, coupling-coefficient, and thermal-bridge metrics;
- analytic and manufactured-solution verification cases;
- plotting and point/profile sampling utilities.

The project is currently an **alpha/research codebase**. SI units are used internally.

## Quick start

The recommended installation uses a Conda environment for the compiled DOLFINx/PETSc/MPI/Gmsh stack and an editable pip install for `toast` itself:

```bash
conda env create -f environment.yml
conda activate toast

# Important on macOS if Homebrew OpenMP paths have been exported globally:
unset DYLD_LIBRARY_PATH

pip install -e .
pytest -q
```

A successful install should currently run the verification and benchmark suite without failures.

## Run examples

Start with the layered wall:

```bash
python examples/layered_wall.py
```

Run the 2-D manufactured-solution example:

```bash
python examples/manufactured.py
```

Run the more complex wall/floor/foundation junction:

```bash
python examples/wall_floor_foundation.py
```

The wall/floor/foundation example uses a continuous exterior-insulation region, a layered above-floor wall, a wood floor/rim region extending into the interior, a concrete foundation, and separate interior/exterior convection boundaries.

## Core architecture

```text
Geometry description
    LayeredWallSpec / RectangleRegion / BoundaryRule
                |
                v
       Gmsh geometry + mesh
                |
                v
   DOLFINx mesh + cell/facet MeshTags
          |                 |
          |                 +----> ConvectionBoundary
          v
   Region tag -> Material mapping
          |
          v
 build_material_fields() -> DG0 k, rho, cp
          |
          v
 solve_steady_heat() / solve_transient_heat()
          |
          v
       ThermalResult
          |
          +----> metrics.py
          +----> postprocess.py
          +----> visualization.py
```

A key design rule is that **geometry identity and material identity are different concepts**. For example, several disconnected geometric rectangles can share one cell tag and map to the same material. Similarly, environmental boundaries are classified independently from material interfaces.

## Package layout

```text
src/toast/
├── materials.py       # Material dataclass and default property database
├── regions.py         # Region and Layer semantic descriptions
├── boundaries.py      # ConvectionBoundary
├── solver.py          # Steady/transient DOLFINx heat solvers
├── metrics.py         # Heat-flow and building-physics metrics
├── postprocess.py     # Sampling and analytic layered-wall utilities
├── verification.py    # Manufactured-solution definitions
├── visualization.py   # PyVista/Matplotlib plotting helpers
├── units.py           # INCH and FOOT helpers
└── mesh/
    ├── layered_wall.py # Explicit conforming layered-wall builder
    ├── patchwork.py    # Generic rectangular patchwork builder
    └── stud_wall.py    # Stud thermal-bridge benchmark geometry
```

## Documentation

More detail is available in [`docs/`](docs/index.md):

- [Installation](docs/installation.md)
- [Architecture](docs/architecture.md)
- [Data structures](docs/data_structures.md)
- [Geometry and meshing](docs/geometry.md)
- [Solvers](docs/solvers.md)
- [Metrics and post-processing](docs/metrics.md)
- [Examples](docs/examples.md)
- [Verification and testing](docs/verification.md)
- [Development notes](docs/development.md)

## Current modeling assumptions and limitations

The current solver assumes isotropic, constant material properties in each tagged region. The patchwork geometry builder is designed primarily for **non-overlapping, axis-aligned rectangular partitions**. Gmsh fragmentation creates conforming interfaces, and generated surfaces are presently associated with regions using center-of-mass classification. That is appropriate for the current wall/stud/foundation examples but is not intended as a general constructive-solid-geometry ownership system.

Thin membranes such as WRBs or air barriers are currently represented either as finite-thickness regions in the layered-wall example or omitted from simplified junction examples. Interface/contact resistances are not yet implemented as first-class boundary/interface models. Likewise, air cavities are presently represented with effective material properties rather than a coupled convection/radiation cavity model.

For 2-D cross-sections, integrated heat flow has units of **W/m**, corresponding to one meter of out-of-plane depth.

## macOS note

DOLFINx, PETSc, MPI, SuperLU_DIST, and OpenMP form a compiled numerical stack. Avoid globally setting Homebrew library search paths while using the Conda environment. In particular, a global value such as

```bash
DYLD_LIBRARY_PATH=/.../homebrew/opt/libomp/lib
```

can cause PETSc imports to fail by loading an incompatible OpenMP runtime. Prefer leaving `DYLD_LIBRARY_PATH` unset inside the `toast` environment.

## Status

The current test suite covers:

- analytic layered thermal resistance and temperature profiles;
- preservation of uniform transient equilibrium;
- steady layered-wall FEM agreement with an analytic solution;
- second-order L2 convergence for a 2-D manufactured solution with P1 elements;
- preservation of that convergence through the Gmsh patchwork/tagging path;
- a heterogeneous stud-wall thermal-bridge benchmark including global energy conservation.
