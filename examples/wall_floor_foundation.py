"""
Wall / floor / foundation junction example.

Simplified 2-D building-envelope thermal bridge with:

    - continuous exterior insulation from wall to foundation
    - layered wall assembly
    - wood floor extending into the interior
    - concrete foundation
    - exterior-temperature space below floor
    - interior-temperature space above floor

Coordinates
-----------
x:
    exterior -> interior

y:
    upward

y = 0:
    top of floor / bottom of wall

This is a steady-state example.
"""

from mpi4py import MPI

from toast.materials import MATERIALS
from toast.regions import Region
from toast.boundaries import ConvectionBoundary

from toast.mesh.patchwork import (
    RectangleRegion,
    BoundaryRule,
    build_rectangular_patchwork_mesh,
)

from toast.solver import solve_steady_heat

from toast.metrics import (
    convection_boundary_heat_flow,
    energy_balance_error,
    thermal_coupling_coefficient,
)

from toast.visualization import (
    plot_material_regions,
    plot_temperature_field,
)


# =====================================================================
# Units
# =====================================================================

INCH = 0.0254


# =====================================================================
# Environmental conditions
# =====================================================================

T_EXTERIOR = 273.15      # K = 0 C
T_INTERIOR = 293.15      # K = 20 C

H_EXTERIOR = 20.0        # W / (m^2 K)
H_INTERIOR = 8.0         # W / (m^2 K)


# =====================================================================
# Physical boundary tags
# =====================================================================

EXTERIOR = 101
INTERIOR = 102


# =====================================================================
# Geometry dimensions
# =====================================================================

# Above-floor wall extent.
#
# This does NOT need to be a full wall height for the junction
# calculation. We just want enough distance for the junction
# disturbance to decay toward 1-D wall behavior.
WALL_HEIGHT = 30.0 * INCH


# Foundation extent below wood floor
FOUNDATION_HEIGHT = 30.0 * INCH


# ---------------------------------------------------------------------
# Wall stack, exterior -> interior
# ---------------------------------------------------------------------

EXT_INS_T = 1.0 * INCH

SHEATHING_T = 0.75 * INCH

CAVITY_T = 5.5 * INCH

SERVICE_T = 1.5 * INCH

DRYWALL_T = 0.625 * INCH


# ---------------------------------------------------------------------
# Floor / foundation
# ---------------------------------------------------------------------

FLOOR_T = 4.0 * INCH

# Distance from exterior insulation face to artificial interior
# truncation of the floor.
FLOOR_RIGHT = 40.0 * INCH

# Concrete stem-wall thickness
FOUNDATION_T = 8.625 * INCH


# ---------------------------------------------------------------------
# Mesh
# ---------------------------------------------------------------------

MESH_SIZE = 0.35 * INCH


# =====================================================================
# Through-wall coordinates
# =====================================================================

x_ext0 = 0.0
x_ext1 = x_ext0 + EXT_INS_T

x_sheathing0 = x_ext1
x_sheathing1 = (
    x_sheathing0
    + SHEATHING_T
)

x_cavity0 = x_sheathing1
x_cavity1 = (
    x_cavity0
    + CAVITY_T
)

x_service0 = x_cavity1
x_service1 = (
    x_service0
    + SERVICE_T
)

x_drywall0 = x_service1
x_drywall1 = (
    x_drywall0
    + DRYWALL_T
)

WALL_RIGHT = x_drywall1


# =====================================================================
# Vertical coordinates
# =====================================================================

y_floor_top = 0.0

y_floor_bottom = (
    y_floor_top
    - FLOOR_T
)

y_foundation_top = y_floor_bottom

y_foundation_bottom = (
    y_foundation_top
    - FOUNDATION_HEIGHT
)

y_wall_bottom = 0.0

y_wall_top = WALL_HEIGHT


# =====================================================================
# Region tags
# =====================================================================

TAG_EXT_INS = 1
TAG_SHEATHING = 2
TAG_CAVITY = 3
TAG_SERVICE = 4
TAG_DRYWALL = 5

TAG_WOOD_FLOOR = 6
TAG_CONCRETE = 7


# =====================================================================
# Geometric regions
#
# IMPORTANT:
#
# These are geometric regions, not material definitions.
#
# The exterior insulation is ONE continuous rectangle extending from
# the bottom of the foundation to the top of the wall.
# =====================================================================

GEOMETRY_REGIONS = [

    # -----------------------------------------------------------------
    # Continuous exterior insulation
    # -----------------------------------------------------------------

    RectangleRegion(
        name="exterior_insulation",
        tag=TAG_EXT_INS,

        x=x_ext0,
        y=y_foundation_bottom,

        width=EXT_INS_T,

        height=(
            y_wall_top
            - y_foundation_bottom
        ),
    ),

    # -----------------------------------------------------------------
    # Above-floor wall stack
    # -----------------------------------------------------------------

    RectangleRegion(
        name="sheathing",
        tag=TAG_SHEATHING,

        x=x_sheathing0,
        y=y_wall_bottom,

        width=SHEATHING_T,
        height=WALL_HEIGHT,
    ),

    RectangleRegion(
        name="cavity_insulation",
        tag=TAG_CAVITY,

        x=x_cavity0,
        y=y_wall_bottom,

        width=CAVITY_T,
        height=WALL_HEIGHT,
    ),

    RectangleRegion(
        name="service_chase",
        tag=TAG_SERVICE,

        x=x_service0,
        y=y_wall_bottom,

        width=SERVICE_T,
        height=WALL_HEIGHT,
    ),

    RectangleRegion(
        name="drywall",
        tag=TAG_DRYWALL,

        x=x_drywall0,
        y=y_wall_bottom,

        width=DRYWALL_T,
        height=WALL_HEIGHT,
    ),

    # -----------------------------------------------------------------
    # Wood structural floor / rim region
    #
    # Starts behind the exterior insulation and extends well into
    # the interior.
    # -----------------------------------------------------------------

    RectangleRegion(
        name="wood_floor",
        tag=TAG_WOOD_FLOOR,

        x=x_ext1,
        y=y_floor_bottom,

        width=(
            FLOOR_RIGHT
            - x_ext1
        ),

        height=FLOOR_T,
    ),

    # -----------------------------------------------------------------
    # Concrete foundation
    #
    # Concrete sits immediately behind the continuous exterior
    # insulation and beneath the wood floor.
    # -----------------------------------------------------------------

    RectangleRegion(
        name="concrete_foundation",
        tag=TAG_CONCRETE,

        x=x_ext1,
        y=y_foundation_bottom,

        width=FOUNDATION_T,
        height=FOUNDATION_HEIGHT,
    ),
]


# =====================================================================
# Thermal region -> material mapping
#
# These names must match keys in MATERIALS.
#
# Based on your current material database, use the distinct exterior
# and interior mineral wool definitions rather than the old generic
# "mineral_wool" name.
# =====================================================================

THERMAL_REGIONS = [

    Region(
        name="exterior_insulation",
        material="mineral_wool_ext",
        tag=TAG_EXT_INS,
    ),

    Region(
        name="sheathing",
        material="plywood",
        tag=TAG_SHEATHING,
    ),

    Region(
        name="cavity_insulation",
        material="mineral_wool_int",
        tag=TAG_CAVITY,
    ),

    Region(
        name="service_chase",
        material="air",
        tag=TAG_SERVICE,
    ),

    Region(
        name="drywall",
        material="drywall",
        tag=TAG_DRYWALL,
    ),

    Region(
        name="wood_floor",
        material="wood",
        tag=TAG_WOOD_FLOOR,
    ),

    Region(
        name="concrete_foundation",
        material="concrete",
        tag=TAG_CONCRETE,
    ),
]


# =====================================================================
# Boundary classification
#
# BoundaryRule is applied ONLY to exposed curves returned by Gmsh,
# not material interfaces.
# =====================================================================

# Predicate tolerance. Geometry is O(1 m), so this is intentionally
# much larger than machine precision and much smaller than any
# construction dimension.
TOL = 1.0e-6


def is_exterior(x, y):
    """
    Exposed surfaces seeing exterior-temperature air.

    Includes:

        1. full left face of continuous exterior insulation
        2. underside of floor beyond the concrete foundation
        3. exposed inside face of the concrete below the floor

    The latter two represent the current assumption that the
    below-floor space is at exterior temperature.
    """

    # -------------------------------------------------------------
    # Exterior face of insulation
    # -------------------------------------------------------------

    exterior_insulation_face = (
        abs(x - x_ext0)
        < TOL
    )

    # -------------------------------------------------------------
    # Underside of floor beyond foundation
    # -------------------------------------------------------------

    foundation_right = (
        x_ext1
        + FOUNDATION_T
    )

    exposed_floor_bottom = (
        abs(
            y - y_floor_bottom
        )
        < TOL

        and

        x
        >= foundation_right - TOL
    )

    # -------------------------------------------------------------
    # Interior-facing side of stem wall below floor
    # -------------------------------------------------------------

    exposed_foundation_face = (
        abs(
            x - foundation_right
        )
        < TOL

        and

        y
        <= y_floor_bottom + TOL
    )

    return (
        exterior_insulation_face
        or exposed_floor_bottom
        or exposed_foundation_face
    )


def is_interior(x, y):
    """
    Exposed surfaces seeing conditioned interior air.

    Includes:

        1. interior face of drywall
        2. top surface of floor extending into the interior
    """

    # -------------------------------------------------------------
    # Interior wall surface
    # -------------------------------------------------------------

    interior_wall_face = (
        abs(
            x - WALL_RIGHT
        )
        < TOL

        and

        y
        >= y_floor_top - TOL
    )

    # -------------------------------------------------------------
    # Interior floor surface
    #
    # Only the portion extending inward from the wall is exposed
    # to room air. Beneath the wall itself, the floor interface is
    # internal to the material domain.
    # -------------------------------------------------------------

    interior_floor_top = (
        abs(
            y - y_floor_top
        )
        < TOL

        and

        x
        >= WALL_RIGHT - TOL
    )

    return (
        interior_wall_face
        or interior_floor_top
    )


BOUNDARY_RULES = [

    BoundaryRule(
        name="exterior",
        tag=EXTERIOR,
        predicate=is_exterior,
    ),

    BoundaryRule(
        name="interior",
        tag=INTERIOR,
        predicate=is_interior,
    ),
]


# =====================================================================
# Run model
# =====================================================================

def main():

    comm = MPI.COMM_WORLD

    # -----------------------------------------------------------------
    # Mesh
    # -----------------------------------------------------------------

    mesh_data = (
        build_rectangular_patchwork_mesh(
            GEOMETRY_REGIONS,
            mesh_size=MESH_SIZE,
            boundary_rules=BOUNDARY_RULES,
            comm=comm,
            model_name=(
                "wall_floor_foundation"
            ),
            write_msh=(
                "wall_floor_foundation.msh"
            ),
        )
    )

    domain = mesh_data.mesh
    cell_tags = mesh_data.cell_tags
    facet_tags = mesh_data.facet_tags

    if comm.rank == 0:

        n_cells = (
            domain.topology.index_map(
                domain.topology.dim
            ).size_global
        )

        print()
        print(
            "Wall/floor/foundation mesh"
        )
        print(
            "=========================="
        )
        print(
            f"Cells: {n_cells}"
        )

        print(
            f"Wall width: "
            f"{WALL_RIGHT / INCH:.3f} in"
        )

        print(
            f"Foundation width: "
            f"{FOUNDATION_T / INCH:.3f} in"
        )

    # -----------------------------------------------------------------
    # Environmental boundary conditions
    # -----------------------------------------------------------------

    boundaries = [

        ConvectionBoundary(
            tag=EXTERIOR,
            h=H_EXTERIOR,
            temperature=T_EXTERIOR,
        ),

        ConvectionBoundary(
            tag=INTERIOR,
            h=H_INTERIOR,
            temperature=T_INTERIOR,
        ),
    ]

    # -----------------------------------------------------------------
    # Steady solve
    # -----------------------------------------------------------------

    result = solve_steady_heat(
        domain,
        cell_tags,
        facet_tags,
        THERMAL_REGIONS,
        MATERIALS,
        boundaries,
    )

    # -----------------------------------------------------------------
    # Global heat-flow metrics
    #
    # Use the convection representation for global conservation.
    # -----------------------------------------------------------------

    q_exterior = (
        convection_boundary_heat_flow(
            domain,
            result.temperature,
            facet_tags,
            EXTERIOR,
            h=H_EXTERIOR,
            ambient_temperature=(
                T_EXTERIOR
            ),
        )
    )

    q_interior = (
        convection_boundary_heat_flow(
            domain,
            result.temperature,
            facet_tags,
            INTERIOR,
            h=H_INTERIOR,
            ambient_temperature=(
                T_INTERIOR
            ),
        )
    )

    balance = energy_balance_error(
        q_exterior,
        q_interior,
    )

    delta_T = (
        T_INTERIOR
        - T_EXTERIOR
    )

    coupling = (
        thermal_coupling_coefficient(
            q_interior,
            delta_T,
        )
    )

    if comm.rank == 0:

        print()
        print(
            "Steady thermal result"
        )
        print(
            "====================="
        )

        print(
            f"Exterior Q: "
            f"{q_exterior: .6f} W/m"
        )

        print(
            f"Interior Q: "
            f"{q_interior: .6f} W/m"
        )

        print(
            f"Energy imbalance: "
            f"{balance:.3e}"
        )

        print(
            f"2-D coupling coefficient: "
            f"{coupling:.6f} W/(m K)"
        )

    # -----------------------------------------------------------------
    # Visualization
    #
    # Keep these after all numerical work because PyVista/VTK has
    # previously been the fragile part of the macOS stack.
    # -----------------------------------------------------------------

    if comm.size == 1:

        plot_material_regions(
            domain,
            cell_tags,
            title=(
                "Wall / floor / foundation "
                "material regions"
            ),
        )

        plot_temperature_field(
            domain,
            result.temperature,
            title=(
                "Wall / floor / foundation "
                "temperature field"
            ),
        )


if __name__ == "__main__":
    main()

