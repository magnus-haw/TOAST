"""
Generic rectangular patchwork mesh builder for DOLFINx.

The domain is assembled from axis-aligned rectangular regions. Gmsh
fragments the rectangles to produce a conforming mesh and physical
groups are assigned from the resulting surfaces.

This is intended for wall/floor/foundation cross sections assembled
from rectangular construction regions.
"""

from dataclasses import dataclass
from collections import defaultdict

import gmsh

from mpi4py import MPI
from dolfinx.io import gmsh as gmshio
from typing import Callable


@dataclass(frozen=True)
class RectangleRegion:
    """
    Axis-aligned rectangular geometric region.

    Parameters
    ----------
    name
        Human-readable name.

    tag
        Physical-region tag used in the DOLFINx cell MeshTags.

        Multiple disconnected rectangles may deliberately share the
        same tag. This is useful when several geometric pieces contain
        the same FEM material region.

    x, y
        Lower-left coordinate [m].

    width, height
        Rectangle dimensions [m].
    """

    name: str
    tag: int

    x: float
    y: float

    width: float
    height: float

@dataclass(frozen=True)
class BoundaryRule:
    """
    Classification rule for an exposed geometric boundary.

    predicate receives the center point (x, y) of an exterior
    Gmsh curve and returns True if that curve belongs to this
    boundary.
    """

    name: str
    tag: int
    predicate: Callable[[float, float], bool]

def _validate_regions(regions, mesh_size):

    regions = list(regions)

    if not regions:
        raise ValueError(
            "At least one RectangleRegion is required."
        )

    if mesh_size <= 0:
        raise ValueError(
            "mesh_size must be positive."
        )

    for region in regions:

        if region.width <= 0:
            raise ValueError(
                f"Region '{region.name}' "
                "has non-positive width."
            )

        if region.height <= 0:
            raise ValueError(
                f"Region '{region.name}' "
                "has non-positive height."
            )

    return regions

def _domain_bounds(regions):

    xmin = min(
        region.x
        for region in regions
    )

    xmax = max(
        region.x + region.width
        for region in regions
    )

    ymin = min(
        region.y
        for region in regions
    )

    ymax = max(
        region.y + region.height
        for region in regions
    )

    return xmin, xmax, ymin, ymax

def _build_fragmented_geometry(
    occ,
    regions,
):

    original_entities = []

    for region in regions:

        surface = occ.addRectangle(
            region.x,
            region.y,
            0.0,
            region.width,
            region.height,
        )

        original_entities.append(
            (2, surface)
        )

    if len(original_entities) > 1:

        occ.fragment(
            [original_entities[0]],
            original_entities[1:],
            removeObject=True,
            removeTool=True,
        )

    occ.synchronize()

def _assign_region_physical_groups(
    occ,
    regions,
    geometry_tol,
):

    resulting_surfaces = (
        gmsh.model.getEntities(dim=2)
    )

    surfaces_by_tag = defaultdict(list)
    names_by_tag = {}
    assigned_surfaces = set()

    for dim, surface_tag in resulting_surfaces:

        cx, cy, _ = occ.getCenterOfMass(
            dim,
            surface_tag,
        )

        matches = []

        for region in regions:

            inside = (
                region.x - geometry_tol
                <= cx
                <= region.x
                + region.width
                + geometry_tol

                and

                region.y - geometry_tol
                <= cy
                <= region.y
                + region.height
                + geometry_tol
            )

            if inside:
                matches.append(region)

        if not matches:

            raise RuntimeError(
                f"Surface {surface_tag} at "
                f"({cx:.6g}, {cy:.6g}) "
                "does not belong to any region."
            )

        physical_tags = {
            region.tag
            for region in matches
        }

        if len(physical_tags) != 1:

            match_text = ", ".join(
                f"{r.name}(tag={r.tag})"
                for r in matches
            )

            raise RuntimeError(
                f"Ambiguous ownership of surface "
                f"{surface_tag}: {match_text}"
            )

        region = matches[0]

        surfaces_by_tag[
            region.tag
        ].append(surface_tag)

        names_by_tag.setdefault(
            region.tag,
            region.name,
        )

        assigned_surfaces.add(
            surface_tag
        )

    expected = {
        tag
        for _, tag in resulting_surfaces
    }

    if assigned_surfaces != expected:

        raise RuntimeError(
            "Some generated surfaces were not "
            "assigned to a region."
        )

    for tag, surfaces in surfaces_by_tag.items():

        gmsh.model.addPhysicalGroup(
            2,
            surfaces,
            tag,
        )

        gmsh.model.setPhysicalName(
            2,
            tag,
            names_by_tag[tag],
        )

def _get_external_boundary_curves():

    surfaces = gmsh.model.getEntities(
        dim=2
    )

    boundary_dimtags = (
        gmsh.model.getBoundary(
            surfaces,
            combined=True,
            oriented=False,
            recursive=False,
        )
    )

    return sorted(
        {
            tag
            for dim, tag in boundary_dimtags
            if dim == 1
        }
    )

def _classify_boundary_curves(
    occ,
    boundary_curves,
    boundary_rules,
):

    curves_by_tag = defaultdict(list)

    rule_names = {
        rule.tag: rule.name
        for rule in boundary_rules
    }

    unclassified = []

    for curve_tag in boundary_curves:

        cx, cy, _ = occ.getCenterOfMass(
            1,
            curve_tag,
        )

        matches = [
            rule
            for rule in boundary_rules
            if rule.predicate(cx, cy)
        ]

        if len(matches) > 1:

            raise RuntimeError(
                f"Boundary curve {curve_tag} "
                f"at ({cx:.6g}, {cy:.6g}) "
                "matched multiple rules: "
                + ", ".join(
                    rule.name
                    for rule in matches
                )
            )

        if not matches:

            unclassified.append(
                curve_tag
            )
            continue

        rule = matches[0]

        curves_by_tag[
            rule.tag
        ].append(curve_tag)

    return (
        curves_by_tag,
        rule_names,
        unclassified,
    )

def _add_boundary_physical_groups(
    curves_by_tag,
    boundary_rules,
):

    for rule in boundary_rules:

        curves = curves_by_tag.get(
            rule.tag,
            [],
        )

        if not curves:

            raise RuntimeError(
                f"No external curves matched "
                f"boundary rule '{rule.name}'."
            )

        gmsh.model.addPhysicalGroup(
            1,
            curves,
            rule.tag,
        )

        gmsh.model.setPhysicalName(
            1,
            rule.tag,
            rule.name,
        )

def _generate_mesh(
    mesh_size,
    write_msh=None,
):

    gmsh.option.setNumber(
        "Mesh.MeshSizeMin",
        mesh_size,
    )

    gmsh.option.setNumber(
        "Mesh.MeshSizeMax",
        mesh_size,
    )

    gmsh.option.setNumber(
        "Mesh.ElementOrder",
        1,
    )

    gmsh.option.setNumber(
        "Mesh.Algorithm",
        6,
    )

    gmsh.model.mesh.generate(2)

    if write_msh is not None:
        gmsh.write(str(write_msh))

def _validate_surface_physical_groups():
    """
    Verify that every 2-D Gmsh surface belongs to exactly one
    2-D physical group.

    DOLFINx requires each top-dimensional cell to have a unique
    physical-region assignment.
    """

    physical_ownership = {}

    for dim, physical_tag in gmsh.model.getPhysicalGroups(dim=2):

        entities = gmsh.model.getEntitiesForPhysicalGroup(
            dim,
            physical_tag,
        )

        for entity in entities:

            if entity in physical_ownership:

                raise RuntimeError(
                    f"Gmsh surface {entity} belongs to more than "
                    f"one 2-D physical group: "
                    f"{physical_ownership[entity]} and "
                    f"{physical_tag}."
                )

            physical_ownership[entity] = physical_tag

    all_surface_ids = {
        tag
        for _, tag in gmsh.model.getEntities(dim=2)
    }

    tagged_surface_ids = set(
        physical_ownership
    )

    if tagged_surface_ids != all_surface_ids:

        missing = (
            all_surface_ids
            - tagged_surface_ids
        )

        extra = (
            tagged_surface_ids
            - all_surface_ids
        )

        message = (
            "Mismatch between Gmsh surfaces and 2-D "
            "physical-group assignments."
        )

        if missing:
            message += (
                f"\nUntagged surfaces: {sorted(missing)}"
            )

        if extra:
            message += (
                f"\nUnknown tagged surfaces: {sorted(extra)}"
            )

        raise RuntimeError(message)

def build_rectangular_patchwork_mesh(
    regions,
    *,
    mesh_size,
    boundary_rules=None,
    comm=MPI.COMM_WORLD,
    model_name="patchwork",
    write_msh=None,
):

    regions = _validate_regions(
        regions,
        mesh_size,
    )

    boundary_rules = list(
        boundary_rules or []
    )

    xmin, xmax, ymin, ymax = (
        _domain_bounds(regions)
    )

    domain_scale = max(
        xmax - xmin,
        ymax - ymin,
        1.0,
    )

    geometry_tol = max(
        1e-10,
        1e-8 * domain_scale,
    )

    rank = 0

    if comm.rank == rank:

        gmsh.initialize()

        try:

            gmsh.model.add(
                model_name
            )

            occ = gmsh.model.occ

            _build_fragmented_geometry(
                occ,
                regions,
            )

            _assign_region_physical_groups(
                occ,
                regions,
                geometry_tol,
            )

            boundary_curves = (
                _get_external_boundary_curves()
            )

            (
                curves_by_tag,
                _,
                unclassified,
            ) = _classify_boundary_curves(
                occ,
                boundary_curves,
                boundary_rules,
            )

            _add_boundary_physical_groups(
                curves_by_tag,
                boundary_rules,
            )

            _validate_surface_physical_groups()

            _generate_mesh(
                mesh_size,
                write_msh,
            )

            mesh_data = (
                gmshio.model_to_mesh(
                    gmsh.model,
                    comm,
                    rank,
                    gdim=2,
                )
            )

        finally:

            gmsh.finalize()

    else:

        gmsh.initialize()

        try:

            mesh_data = (
                gmshio.model_to_mesh(
                    gmsh.model,
                    comm,
                    rank,
                    gdim=2,
                )
            )

        finally:

            gmsh.finalize()

    return mesh_data







