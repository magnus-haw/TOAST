# toast documentation

`toast` provides reusable geometry, solver, verification, and post-processing tools for 2-D building-envelope heat-transfer models using Gmsh and DOLFINx.

The package is intentionally split into layers. Geometry code creates a tagged mesh; semantic regions map cell tags to materials; boundary descriptions map facet tags to thermal environments; solvers assemble and solve the heat equation; metrics and visualization operate on the resulting `ThermalResult`.

## Documentation map

| Page | Purpose |
|---|---|
| [Installation](installation.md) | Create the Conda environment, install the package, run tests, and troubleshoot macOS library issues. |
| [Architecture](architecture.md) | Understand module responsibilities and how data flow through the package. |
| [Data structures](data_structures.md) | Reference the main dataclasses and DOLFINx mesh/tag objects. |
| [Geometry and meshing](geometry.md) | Build layered and rectangular patchwork geometries and classify boundaries. |
| [Solvers](solvers.md) | Governing equations, steady/transient APIs, material fields, boundary conditions, and outputs. |
| [Metrics and post-processing](metrics.md) | Heat flow, energy balance, U-values, coupling coefficients, psi values, sampling, and plotting. |
| [Examples](examples.md) | Run and interpret the supplied examples. |
| [Verification](verification.md) | Numerical verification strategy and current tests. |
| [Development](development.md) | Add materials, geometries, solvers, metrics, and tests without breaking package abstractions. |

## Units and conventions

All geometry and physics use SI units internally:

- length: m
- temperature: K
- conductivity: W/(m K)
- density: kg/m^3
- specific heat: J/(kg K)
- convection coefficient: W/(m^2 K)
- 2-D integrated heat flow: W/m of out-of-plane depth

Convenience constants `toast.units.INCH` and `toast.units.FOOT` are provided for geometry definitions.

The coordinate convention is geometry-specific. The supplied wall examples generally use increasing `x` from exterior to interior. Boundary tags are semantic; do not assume that an `EXTERIOR` boundary always coincides with the global minimum-x boundary in complex junction geometries.
