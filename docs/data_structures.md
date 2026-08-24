# Data structures

This page summarizes the primary data structures passed between package modules.

## `Material`

```python
from toast.materials import Material

wood = Material(
    k=0.12,
    rho=450.0,
    cp=1600.0,
)
```

Fields:

| Field | Meaning | Units |
|---|---|---|
| `k` | isotropic thermal conductivity | W/(m K) |
| `rho` | density | kg/m^3 |
| `cp` | specific heat capacity | J/(kg K) |

The default database is `toast.materials.MATERIALS`. Properties are currently constant; temperature-dependent properties are not yet represented.

## `Region`

```python
from toast.regions import Region

cavity = Region(
    name="cavity_insulation",
    material="mineral_wool_int",
    tag=3,
)
```

`Region` is the semantic mapping between a DOLFINx/Gmsh cell tag and a material name:

```text
cell tag 3 -> Region("cavity_insulation", "mineral_wool_int", 3)
           -> MATERIALS["mineral_wool_int"]
           -> k, rho, cp on those cells
```

The solver uses `Region.tag` to locate cells and `Region.material` to look up properties.

## `Layer`

`Layer` extends `Region`:

```python
from toast.regions import Layer
from toast.units import INCH

layer = Layer(
    name="sheathing",
    material="plywood",
    tag=5,
    thickness=0.75 * INCH,
    n_through=4,
)
```

`n_through` is used by the layered-wall mesher to choose a nominal point mesh size through that layer.

## `RectangleRegion`

```python
from toast.mesh.patchwork import RectangleRegion

foundation = RectangleRegion(
    name="concrete_foundation",
    tag=7,
    x=0.0254,
    y=-0.86,
    width=0.22,
    height=0.76,
)
```

This is a **geometric** description, not a material description. It defines an axis-aligned rectangle and the cell tag that its generated surfaces should receive.

Several disconnected rectangles may intentionally share the same tag. The stud-wall benchmark uses this to assign both upper and lower insulation rectangles to one insulation tag.

## `BoundaryRule`

```python
from toast.mesh.patchwork import BoundaryRule

EXTERIOR = 101

def is_exterior(x, y):
    return abs(x) < 1e-6

rule = BoundaryRule(
    name="exterior",
    tag=EXTERIOR,
    predicate=is_exterior,
)
```

The predicate receives the center-of-mass coordinate `(x, y)` of an **exposed external Gmsh curve**. It returns `True` when the curve should belong to that physical boundary.

A curve must not match more than one boundary rule. Exposed curves that match no rule remain untagged and therefore receive the natural PDE boundary condition.

## `ConvectionBoundary`

```python
from toast.boundaries import ConvectionBoundary

bc = ConvectionBoundary(
    tag=101,
    h=20.0,
    temperature=273.15,
    name="exterior",
)
```

This represents the Robin condition

```text
-k grad(T) . n = h (T - T_inf)
```

with `temperature` equal to `T_inf`.

## `LayeredWallSpec`

```python
from toast.mesh.layered_wall import LayeredWallSpec

spec = LayeredWallSpec(
    layers=my_layers,
    height=0.5,
    vertical_mesh_size=0.02,
)
```

Important attributes include `layers`, `height`, `vertical_mesh_size`, `boundary_tags`, and the derived `width` property.

## `StudWallSpec`

`StudWallSpec` defines the benchmark section dimensions, material/tag identifiers, and mesh size for the stud-wall geometry.

## DOLFINx `MeshData`

Both mesh builders return the object produced by `dolfinx.io.gmsh.model_to_mesh`. The package uses:

```python
mesh_data.mesh
mesh_data.cell_tags
mesh_data.facet_tags
mesh_data.physical_groups
```

`cell_tags` identify material regions. `facet_tags` identify environmental boundaries.

## `SolverConfig`

For transient solves:

```python
from toast.solver import SolverConfig

config = SolverConfig(
    initial_temperature=293.15,
    dt=60.0,
    t_end=24 * 3600.0,
    volumetric_source=0.0,
    degree=1,
    output_interval=600.0,
    output_path="wall_temperature.bp",
    petsc_options_prefix="thermal_",
)
```

The current transient solver assumes a scalar, spatially uniform initial temperature and scalar volumetric source.

## `ThermalResult`

Both solver paths return:

```python
@dataclass
class ThermalResult:
    temperature
    function_space
    conductivity
    density
    heat_capacity
    time
```

The property fields are DG0 functions. Temperature uses a continuous Lagrange function space of the requested degree.
