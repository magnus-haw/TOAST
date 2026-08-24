# Heat solvers

## Governing equation

The transient model solves heterogeneous heat conduction of the form

```text
rho cp dT/dt - div(k grad(T)) = Q
```

where `k`, `rho`, and `cp` are piecewise-constant DG0 fields assembled from cell tags and the material database.

The steady model removes the storage term:

```text
-div(k grad(T)) = Q
```

The current implementation assumes isotropic scalar conductivity.

## Convection boundaries

`ConvectionBoundary` represents

```text
-k grad(T) . n = h (T - T_inf)
```

The Robin term enters the weak form as a boundary contribution to both the bilinear and linear forms.

Boundaries that are not explicitly included in the solver are left as the natural zero-flux Neumann condition.

## Material-field assembly

`build_material_fields(domain, cell_tags, regions, materials)` creates DG0 functions for:

- thermal conductivity `k`;
- density `rho`;
- specific heat `cp`.

For every semantic `Region`, the solver finds cells with `cell_tags.find(region.tag)` and fills those property arrays from `materials[region.material]`.

The function initializes all values to NaN and raises an error if any cell remains unassigned. It also raises a `KeyError` when a region references an unknown material name.

## Steady solve

```python
from toast.solver import solve_steady_heat

result = solve_steady_heat(
    domain,
    cell_tags,
    facet_tags,
    regions,
    materials,
    boundaries,
    volumetric_source=0.0,
    degree=1,
    petsc_options_prefix="thermal_steady_",
)
```

The current PETSc configuration is:

```text
KSP: CG
PC:  GAMG
relative tolerance: 1e-10
absolute tolerance: 1e-12
error if not converged: true
```

The heat-conduction matrix with positive conductivity and convection terms is expected to be symmetric positive definite under the supported configurations, making CG a reasonable default.

## Transient solve

```python
from toast.solver import SolverConfig, solve_transient_heat

config = SolverConfig(
    initial_temperature=293.15,
    dt=60.0,
    t_end=24 * 3600.0,
    output_interval=600.0,
    output_path="wall_temperature.bp",
)

result = solve_transient_heat(
    domain,
    cell_tags,
    facet_tags,
    regions,
    materials,
    boundaries,
    config=config,
)
```

Time integration uses backward Euler:

```text
rho cp (T[n+1] - T[n]) / dt
    - div(k grad(T[n+1])) = Q
```

The solver writes VTX/ADIOS2 output when `output_path` is not `None`. `output_interval=None` disables periodic output scheduling, and `output_path=None` disables the writer entirely.

## `ThermalResult`

Both solvers return a `ThermalResult` containing:

```python
result.temperature
result.function_space
result.conductivity
result.density
result.heat_capacity
result.time
```

For a steady solve, `time` is `0.0`. For a transient solve, it is the final simulation time reached.

## Solver independence from geometry

The solver does not inspect geometric coordinates to decide materials or boundary identities. It requires only:

- a DOLFINx mesh;
- cell tags;
- facet tags;
- semantic region-to-material mappings;
- boundary-condition objects.

This is what allows the same solver to operate on layered walls, stud walls, and wall/floor/foundation junctions.

## Current limitations

The current solver API supports:

- constant scalar material properties;
- scalar uniform volumetric source values;
- constant convection coefficients and ambient temperatures per tagged boundary;
- scalar uniform initial temperature in transient runs.

Temperature-dependent materials, anisotropic conductivity, contact/interface resistance, radiation, nonlinear boundary conditions, and spatial/time-dependent sources are not yet first-class APIs.
