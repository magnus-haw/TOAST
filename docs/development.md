# Development guide

## Repository workflow

Activate the Conda environment and install the repository in editable mode:

```bash
conda activate toast
unset DYLD_LIBRARY_PATH  # macOS, if needed
pip install -e .
```

Run tests after changes:

```bash
pytest -q
```

## Adding a material

Add a `Material` entry in `src/toast/materials.py`:

```python
MATERIALS["my_material"] = Material(
    k=0.15,
    rho=500.0,
    cp=1200.0,
)
```

Keep material names descriptive of the physical material, not its geometric role. For example, prefer `wood` over `rim_joist` if the same property set is used for studs, plates, and floor framing.

If two product classes genuinely have different properties, use separate keys such as the existing `mineral_wool_ext` and `mineral_wool_int`.

## Adding a semantic region

A geometric cell tag must be mapped to a material before solving:

```python
Region(
    name="rim_region",
    material="wood",
    tag=12,
)
```

Every cell must receive exactly one semantic material mapping. `build_material_fields()` will fail if cells remain NaN.

## Adding a rectangular patchwork geometry

1. Define non-overlapping `RectangleRegion` objects.
2. Reuse a physical tag for disconnected geometric pieces that should share one semantic region.
3. Define `BoundaryRule` predicates only for environmental boundaries needed by the PDE.
4. Build with `build_rectangular_patchwork_mesh()`.
5. Plot material tags before trusting thermal results.
6. Check energy conservation.
7. Perform mesh/domain-size convergence for production metrics.

Example:

```python
mesh_data = build_rectangular_patchwork_mesh(
    geometry_regions,
    mesh_size=0.01,
    boundary_rules=boundary_rules,
    model_name="my_junction",
)
```

## Adding a boundary condition type

At present the solver knows only `ConvectionBoundary`. A new boundary condition should generally have:

- a small semantic dataclass in `boundaries.py`;
- assembly logic in `solver.py`;
- a corresponding integrated-flux utility in `metrics.py` when global conservation requires it;
- an isolated verification test.

Avoid embedding geometric coordinate logic inside the solver. Geometry should provide facet tags; physics should consume them.

## Adding solver physics

Changes such as temperature-dependent conductivity, anisotropic conductivity, radiation, contact resistance, or nonlinear convection should be introduced behind the geometry-independent solver API where possible.

Before extending the production geometry, first create a small verification problem that exercises the new mathematical term with a known or strongly constrained solution.

## Adding metrics

Global building-performance quantities belong in `metrics.py`. Point/profile operations and analytic solution helpers belong in `postprocess.py`.

There is currently duplicated legacy functionality between these two modules for some quantities such as boundary heat flow and U-value helpers. New code should prefer `metrics.py`; the duplicate functions in `postprocess.py` can be deprecated/removed in a later cleanup once callers are migrated.

## Gmsh physical groups

DOLFINx expects every top-dimensional mesh cell to be represented by one physical group. Keep the defensive `_validate_surface_physical_groups()` check in the patchwork builder.

For the patchwork builder, environmental boundary classification should only operate on the combined external Gmsh boundary. Do not classify all Gmsh curves directly, because internal material interfaces would then be candidates for environmental BCs.

## Current geometry boundary of applicability

The current `RectangleRegion`/centroid-ownership architecture is good for assembly partitions made from rectangles. Do not stretch it into a general overlapping CAD system.

When requirements include embedded soil/foundation subtraction, curved/nonrectangular regions, holes, or competing overlapping material primitives, refactor semantic ownership to use Gmsh's Boolean-fragment output mapping rather than centroid containment.

## Code quality

The `pyproject.toml` includes configuration for pytest, Ruff, mypy, and coverage. The scientific libraries do not provide complete type information, so mypy is deliberately permissive around external imports.
