# Architecture

## Design goals

The package separates **geometry**, **material semantics**, **boundary semantics**, **numerical solution**, and **analysis**. This is important for building-envelope models because the same physical material may occur in many disconnected geometric locations, and the same environmental boundary may consist of several disconnected exterior curves.

The intended data flow is:

```text
Geometry specification
   |-- LayeredWallSpec / Layer
   |-- RectangleRegion / BoundaryRule
   v
Gmsh CAD + physical groups
   v
DOLFINx MeshData
   |-- mesh
   |-- cell_tags
   |-- facet_tags
   v
Semantic Region definitions + Material database
   v
DG0 property fields: k, rho, cp
   v
ConvectionBoundary objects
   v
steady/transient solver
   v
ThermalResult
   |-- temperature
   |-- function_space
   |-- conductivity
   |-- density
   |-- heat_capacity
   `-- time
   v
metrics / post-processing / visualization
```

## Module responsibilities

### `materials.py`

Defines `Material(k, rho, cp)` and the default `MATERIALS` mapping. Material names are semantic identifiers such as `wood`, `concrete`, `mineral_wool_ext`, and `mineral_wool_int`.

### `regions.py`

Defines `Region(name, material, tag)`, which links a cell tag to a material name, and `Layer`, which extends `Region` with thickness and mesh-resolution metadata for layered walls.

A region is not a CAD primitive. It is the semantic interpretation of a cell tag.

### `boundaries.py`

Defines `ConvectionBoundary(tag, h, temperature, name="")`. The solver uses the tag to integrate the Robin boundary term over the corresponding facet physical group.

### `mesh/layered_wall.py`

Builds a special-purpose layered wall using shared Gmsh line entities at all material interfaces. This is a useful robust path for nominally 1-D walls and analytic verification.

### `mesh/patchwork.py`

Builds a more general non-overlapping rectangular assembly. `RectangleRegion` objects create CAD rectangles; Gmsh fragmentation makes interfaces conforming; cell physical groups are assigned by region tag; `BoundaryRule` objects classify exposed outer curves.

The public function is `build_rectangular_patchwork_mesh(...)`.

### `mesh/stud_wall.py`

Provides a small benchmark geometry demonstrating disconnected insulation rectangles sharing one material-region tag and a wood stud acting as a thermal bridge.

### `solver.py`

Contains material-field construction and geometry-independent steady/transient heat solvers. The solver knows about tags and materials but does not know whether a mesh came from the layered-wall builder, the patchwork builder, or another future geometry source.

### `metrics.py`

Contains integrated boundary heat-flow calculations and higher-level building-physics quantities such as energy imbalance, U-value references, thermal coupling coefficient, and linear thermal transmittance.

### `postprocess.py`

Contains point sampling and layered-wall analytic helpers. There is currently some functional overlap between `postprocess.py` and `metrics.py`; new global thermal-performance metrics should generally be added to `metrics.py`.

### `verification.py`

Defines the exact temperature and source expressions used by the manufactured-solution tests.

### `visualization.py`

Contains PyVista mesh/field views and Matplotlib temperature-profile helpers. Visualization is intentionally kept outside the solver.

## Important invariants

The meshing code should preserve these invariants:

1. Every top-dimensional cell belongs to exactly one material physical group.
2. Material interfaces are conforming mesh interfaces, not coincident disconnected surfaces.
3. Environmental facet groups contain only exposed outer boundaries, not internal material interfaces.
4. Multiple disconnected CAD regions may share the same physical cell tag if they represent the same semantic material region.
5. Unclassified exposed boundaries are left to the natural FEM boundary condition; for the current heat equation this means zero normal conductive flux.

## Why boundary predicates are separate from region geometry

In a simple wall, exterior and interior can be identified as minimum-x and maximum-x boundaries. A wall/floor/foundation junction breaks that assumption: the exterior boundary can include a vertical insulation face, a floor underside, and a foundation face, while the interior can include both the drywall face and floor top.

`BoundaryRule` therefore classifies each exposed Gmsh curve by a predicate evaluated at the curve center of mass. The geometry builder first asks Gmsh for the **combined external boundary**, so internal material interfaces are excluded before classification.
