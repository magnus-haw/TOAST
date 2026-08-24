# src/toast/mesh/stud_wall.py

from dataclasses import dataclass

from mpi4py import MPI

from .patchwork import (
    RectangleRegion,
    BoundaryRule,
    build_rectangular_patchwork_mesh,
)


INCH = 0.0254


EXTERIOR = 101
INTERIOR = 102
BOTTOM = 103
TOP = 104


@dataclass(frozen=True)
class StudWallSpec:
    thickness: float = 5.5 * INCH
    height: float = 16.0 * INCH
    stud_width: float = 1.5 * INCH

    insulation_tag: int = 1
    stud_tag: int = 2

    insulation_material: str = "mineral_wool_int"
    stud_material: str = "wood"

    mesh_size: float = 0.25 * INCH


def stud_wall_regions(
    spec: StudWallSpec,
):
    """
    Return three rectangular regions:

        upper insulation
        wood stud
        lower insulation
    """

    remaining = (
        spec.height
        - spec.stud_width
    )

    lower_height = remaining / 2.0
    upper_height = remaining / 2.0

    return [
        RectangleRegion(
            name="lower_insulation",
            tag=spec.insulation_tag,
            x=0.0,
            y=0.0,
            width=spec.thickness,
            height=lower_height,
        ),

        RectangleRegion(
            name="stud",
            tag=spec.stud_tag,
            x=0.0,
            y=lower_height,
            width=spec.thickness,
            height=spec.stud_width,
        ),

        RectangleRegion(
            name="upper_insulation",
            tag=spec.insulation_tag,
            x=0.0,
            y=lower_height + spec.stud_width,
            width=spec.thickness,
            height=upper_height,
        ),
    ]


def build_stud_wall_mesh(
    spec=StudWallSpec(),
    *,
    comm=MPI.COMM_WORLD,
):
    """
    Build the stud-wall benchmark mesh.

    Boundary convention
    -------------------
    x = 0
        exterior

    x = spec.thickness
        interior

    y = 0
        bottom / adiabatic

    y = spec.height
        top / adiabatic
    """

    regions = stud_wall_regions(
        spec
    )

    # Geometry is O(0.1 m), so this is generous relative to
    # OpenCASCADE geometric tolerances.
    tol = 1e-8

    boundary_rules = [

        BoundaryRule(
            name="exterior",
            tag=EXTERIOR,
            predicate=lambda x, y: (
                abs(x - 0.0) < tol
            ),
        ),

        BoundaryRule(
            name="interior",
            tag=INTERIOR,
            predicate=lambda x, y: (
                abs(
                    x - spec.thickness
                ) < tol
            ),
        ),

    ]

    return build_rectangular_patchwork_mesh(
        regions,
        mesh_size=spec.mesh_size,
        boundary_rules=boundary_rules,
        comm=comm,
        model_name="stud_wall",
    )

