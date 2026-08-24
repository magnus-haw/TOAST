"""Conforming Gmsh mesh for a layered 2-D wall cross-section.

Each material layer is a separate plane surface, adjacent layers share the
same geometric interface line, and each top-dimensional cell belongs to one
and only one physical group. This preserves the topology of the original
working wall_mesh.py while separating geometry from the thermal solver.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import gmsh
from mpi4py import MPI
from dolfinx.io import gmsh as gmshio

from toast.regions import Layer
from toast.units import FOOT, INCH


@dataclass(frozen=True)
class BoundaryTags:
    exterior: int = 101
    interior: int = 102
    bottom: int = 103
    top: int = 104


EXTERIOR = 101
INTERIOR = 102
BOTTOM = 103
TOP = 104


DEFAULT_LAYERS = (
    Layer("exterior_render", "fiber_cement", 1, 3.0 / 8.0 * INCH, 3),
    Layer("drainage_gap", "air", 2, 1.0 / 16.0 * INCH, 2),
    Layer("exterior_insulation", "mineral_wool_ext", 3, 1.0 * INCH, 5),
    Layer("WRB", "air", 4, 1.0 / 16.0 * INCH, 2),
    Layer("sheathing", "plywood", 5, 3.0 / 4.0 * INCH, 4),
    Layer("framing", "mineral_wool_int", 6, 5.5 * INCH, 10),
    Layer("air_barrier", "air", 7, 1.0 / 16.0 * INCH, 2),
    Layer("service_chase", "air", 8, 1.5 * INCH, 5),
    Layer("drywall", "drywall", 9, 0.5 * INCH, 3),
)


@dataclass(frozen=True)
class LayeredWallSpec:
    layers: Sequence[Layer] = DEFAULT_LAYERS
    height: float = 1.0 * FOOT
    vertical_mesh_size: float = 1.0 * INCH
    boundary_tags: BoundaryTags = BoundaryTags()

    @property
    def width(self) -> float:
        return sum(layer.thickness for layer in self.layers)


DEFAULT_WALL = LayeredWallSpec()


def build_layered_wall_mesh(
    spec: LayeredWallSpec = DEFAULT_WALL,
    *,
    comm=MPI.COMM_WORLD,
    msh_path: str | Path | None = None,
    verbose: bool = True,
):
    """Generate a conforming layered-wall Gmsh model and import it to DOLFINx."""

    rank = 0
    layers = tuple(spec.layers)
    tags = spec.boundary_tags

    if comm.rank == rank:
        gmsh.initialize()
        gmsh.model.add("layered_wall")
        geo = gmsh.model.geo

        x_coords = [0.0]
        for layer in layers:
            x_coords.append(x_coords[-1] + layer.thickness)

        bottom_points = []
        top_points = []
        for x in x_coords:
            bottom_points.append(
                geo.addPoint(x, 0.0, 0.0, spec.vertical_mesh_size)
            )
            top_points.append(
                geo.addPoint(x, spec.height, 0.0, spec.vertical_mesh_size)
            )

        # One shared vertical entity at every material interface.
        vertical_lines = [
            geo.addLine(bottom_points[i], top_points[i])
            for i in range(len(x_coords))
        ]

        bottom_lines = []
        top_lines = []
        for i in range(len(layers)):
            bottom_lines.append(
                geo.addLine(bottom_points[i], bottom_points[i + 1])
            )
            top_lines.append(
                geo.addLine(top_points[i + 1], top_points[i])
            )

        surfaces = []
        for i in range(len(layers)):
            loop = geo.addCurveLoop(
                [
                    bottom_lines[i],
                    vertical_lines[i + 1],
                    top_lines[i],
                    -vertical_lines[i],
                ]
            )
            surfaces.append(geo.addPlaneSurface([loop]))

        geo.synchronize()

        for layer, surface in zip(layers, surfaces):
            gmsh.model.addPhysicalGroup(2, [surface], tag=layer.tag)
            gmsh.model.setPhysicalName(2, layer.tag, layer.name)

        gmsh.model.addPhysicalGroup(1, [vertical_lines[0]], tag=tags.exterior)
        gmsh.model.setPhysicalName(1, tags.exterior, "exterior")
        gmsh.model.addPhysicalGroup(1, [vertical_lines[-1]], tag=tags.interior)
        gmsh.model.setPhysicalName(1, tags.interior, "interior")
        gmsh.model.addPhysicalGroup(1, bottom_lines, tag=tags.bottom)
        gmsh.model.setPhysicalName(1, tags.bottom, "bottom")
        gmsh.model.addPhysicalGroup(1, top_lines, tag=tags.top)
        gmsh.model.setPhysicalName(1, tags.top, "top")

        # Preserve the original layer-aware point sizing.
        for i, layer in enumerate(layers):
            dx = layer.thickness / layer.n_through
            gmsh.model.mesh.setSize(
                [
                    (0, bottom_points[i]),
                    (0, top_points[i]),
                    (0, bottom_points[i + 1]),
                    (0, top_points[i + 1]),
                ],
                dx,
            )

        gmsh.option.setNumber("Mesh.ElementOrder", 1)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.model.mesh.generate(2)

        # Defensive check: DOLFINx requires one top-dimensional physical tag/cell.
        seen_surfaces: dict[int, str] = {}
        for dim, physical_tag in gmsh.model.getPhysicalGroups(2):
            name = gmsh.model.getPhysicalName(dim, physical_tag)
            entities = gmsh.model.getEntitiesForPhysicalGroup(dim, physical_tag)
            for entity in entities:
                if int(entity) in seen_surfaces:
                    raise RuntimeError(
                        f"Surface {entity} is tagged by both "
                        f"{seen_surfaces[int(entity)]!r} and {name!r}."
                    )
                seen_surfaces[int(entity)] = name

        if len(seen_surfaces) != len(layers):
            raise RuntimeError(
                f"Expected {len(layers)} tagged surfaces, "
                f"found {len(seen_surfaces)}."
            )

        if msh_path is not None:
            gmsh.write(str(msh_path))

        if verbose:
            print(f"Total wall thickness: {spec.width / INCH:.4f} in")
            print(f"Wall height: {spec.height / FOOT:.3f} ft")

    mesh_data = gmshio.model_to_mesh(gmsh.model, comm, rank, gdim=2)

    if comm.rank == rank:
        gmsh.finalize()

    return mesh_data


# Backward-friendly alias for the original script naming.
build_wall_mesh = build_layered_wall_mesh
