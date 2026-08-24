# Geometry and meshing

## Geometry backends

`toast` currently exposes two meshing strategies:

- `build_layered_wall_mesh()` for ordered planar wall layers;
- `build_rectangular_patchwork_mesh()` for non-overlapping axis-aligned rectangular regions.

Both return DOLFINx `MeshData` imported from Gmsh physical groups.

## Layered walls

A layered wall is described by `Layer` objects inside a `LayeredWallSpec`:

```python
from toast.mesh.layered_wall import LayeredWallSpec, build_layered_wall_mesh
from toast.regions import Layer
from toast.units import INCH

layers = (
    Layer("exterior_insulation", "mineral_wool_ext", 1, 1.0 * INCH, 5),
    Layer("sheathing", "plywood", 2, 0.75 * INCH, 4),
    Layer("cavity", "mineral_wool_int", 3, 5.5 * INCH, 10),
    Layer("drywall", "drywall", 4, 0.5 * INCH, 3),
)

spec = LayeredWallSpec(
    layers=layers,
    height=1.0,
    vertical_mesh_size=0.025,
)

mesh_data = build_layered_wall_mesh(spec)
```

The builder explicitly creates one shared vertical Gmsh line at each material interface. This avoids coincident disconnected boundaries and guarantees a conforming mesh between adjacent layers.

The standard boundary-tag convention is stored in `BoundaryTags` and defaults to exterior=101, interior=102, bottom=103, top=104.

## Rectangular patchworks

A patchwork is a collection of `RectangleRegion` objects:

```python
from toast.mesh.patchwork import RectangleRegion

regions = [
    RectangleRegion("insulation", 1, x=0.0, y=0.0, width=0.14, height=0.30),
    RectangleRegion("wood_floor", 2, x=0.14, y=-0.10, width=0.8, height=0.10),
]
```

The public builder is:

```python
mesh_data = build_rectangular_patchwork_mesh(
    regions,
    mesh_size=0.01,
    boundary_rules=boundary_rules,
    model_name="assembly",
    write_msh="assembly.msh",
)
```

Internally it performs the following steps:

1. validate region dimensions and mesh size;
2. create Gmsh OpenCASCADE rectangles;
3. fragment the rectangles together so shared interfaces are conforming;
4. associate generated surfaces with semantic region tags;
5. add one 2-D physical group per unique cell tag;
6. obtain the combined exposed boundary from Gmsh;
7. classify exposed curves with `BoundaryRule` predicates;
8. create facet physical groups;
9. verify that every generated surface belongs to exactly one 2-D physical group;
10. generate the triangular mesh and import it into DOLFINx.

## Boundary rules

For complex junctions, exterior/interior boundaries are not simply domain extrema. A rule can describe disconnected boundary pieces:

```python
TOL = 1e-6
EXTERIOR = 101

foundation_right = 0.25
floor_bottom = -0.10

def is_exterior(x, y):
    left_face = abs(x) < TOL
    floor_underside = abs(y - floor_bottom) < TOL and x >= foundation_right - TOL
    foundation_face = abs(x - foundation_right) < TOL and y <= floor_bottom + TOL
    return left_face or floor_underside or foundation_face
```

The patchwork builder first asks Gmsh for the combined outer boundary. Internal material interfaces therefore never reach the boundary-rule classifier.

Unclassified exposed boundaries are allowed. With the current heat equation, omitting a boundary term gives the natural zero-flux Neumann condition. This is useful for artificial truncation/symmetry boundaries.

## Material tags versus geometry names

`RectangleRegion.tag` is the key connection between geometry and thermal semantics. Geometry names are for readability; the solver uses integer tags.

For example, disconnected rectangles can share a single tag:

```python
RectangleRegion("lower_insulation", tag=1, ...)
RectangleRegion("upper_insulation", tag=1, ...)
```

and one semantic region can map tag 1 to `mineral_wool_int`:

```python
Region("insulation", "mineral_wool_int", tag=1)
```

## Current patchwork limitations

The patchwork builder is designed for a non-overlapping rectangular partition. After fragmentation, it currently classifies each generated surface using its center of mass and the original rectangle extents. That is robust for the existing partition-style geometries but becomes ambiguous for general overlapping constructive-solid-geometry inputs.

Do not use rectangle ordering as an implicit material-precedence mechanism. If future models use soil backgrounds with embedded/overlapping foundation objects, nonrectangular regions, holes, or boolean subtraction, the geometry ownership model should be upgraded to propagate semantic identity through the Gmsh fragment mapping.

## Mesh resolution

The layered-wall builder has layer-specific `n_through` information. The patchwork builder currently uses a single global `mesh_size` through `Mesh.MeshSizeMin` and `Mesh.MeshSizeMax`.

For junction models, convergence should be checked with respect to both:

- mesh resolution;
- artificial domain extent (wall height, foundation depth, and floor extension).

Local refinement around material junctions is a likely future extension.
