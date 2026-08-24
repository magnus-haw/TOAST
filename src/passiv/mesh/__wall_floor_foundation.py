"""
Wall / floor / foundation thermal-bridge example.

Coordinates
-----------
x:
    exterior -> interior

y:
    upward

y = 0:
    top of structural floor / bottom of wall
"""

from toast.mesh.patchwork import RectangleRegion


INCH = 0.0254


# ============================================================
# Dimensions
# ============================================================

WALL_HEIGHT = 30.0 * INCH
FOUNDATION_HEIGHT = 30.0 * INCH

EXTERIOR_INSULATION_T = 1.0 * INCH
SHEATHING_T = 0.75 * INCH
CAVITY_T = 5.5 * INCH
SERVICE_T = 1.5 * INCH
DRYWALL_T = 0.625 * INCH

FLOOR_T = 4.0 * INCH
FLOOR_LENGTH = 30.0 * INCH

FOUNDATION_T = 8.625 * INCH


# ============================================================
# Derived wall coordinates
# ============================================================

x0 = 0.0

x_ext_ins_0 = x0
x_ext_ins_1 = x_ext_ins_0 + EXTERIOR_INSULATION_T

x_sheathing_0 = x_ext_ins_1
x_sheathing_1 = x_sheathing_0 + SHEATHING_T

x_cavity_0 = x_sheathing_1
x_cavity_1 = x_cavity_0 + CAVITY_T

x_service_0 = x_cavity_1
x_service_1 = x_service_0 + SERVICE_T

x_drywall_0 = x_service_1
x_drywall_1 = x_drywall_0 + DRYWALL_T

WALL_WIDTH = x_drywall_1


# ============================================================
# Material tags
# ============================================================

TAG_EXT_INS = 1
TAG_SHEATHING = 2
TAG_CAVITY_INSULATION = 3
TAG_SERVICE_CHASE = 4
TAG_DRYWALL = 5

TAG_WOOD_FLOOR = 6
TAG_CONCRETE = 7


# ============================================================
# Regions
# ============================================================

REGIONS = [

    # --------------------------------------------------------
    # Wall above floor
    # --------------------------------------------------------

    RectangleRegion(
        name="exterior_insulation",
        tag=TAG_EXT_INS,
        x=x_ext_ins_0,
        y=0.0,
        width=EXTERIOR_INSULATION_T,
        height=WALL_HEIGHT,
    ),

    RectangleRegion(
        name="sheathing",
        tag=TAG_SHEATHING,
        x=x_sheathing_0,
        y=0.0,
        width=SHEATHING_T,
        height=WALL_HEIGHT,
    ),

    RectangleRegion(
        name="cavity_insulation",
        tag=TAG_CAVITY_INSULATION,
        x=x_cavity_0,
        y=0.0,
        width=CAVITY_T,
        height=WALL_HEIGHT,
    ),

    RectangleRegion(
        name="service_chase",
        tag=TAG_SERVICE_CHASE,
        x=x_service_0,
        y=0.0,
        width=SERVICE_T,
        height=WALL_HEIGHT,
    ),

    RectangleRegion(
        name="drywall",
        tag=TAG_DRYWALL,
        x=x_drywall_0,
        y=0.0,
        width=DRYWALL_T,
        height=WALL_HEIGHT,
    ),

    # --------------------------------------------------------
    # Wood floor / rim region
    #
    # Begins behind the exterior insulation and extends
    # significantly into the building interior.
    # --------------------------------------------------------

    RectangleRegion(
        name="wood_floor",
        tag=TAG_WOOD_FLOOR,
        x=x_ext_ins_1,
        y=-FLOOR_T,
        width=FLOOR_LENGTH,
        height=FLOOR_T,
    ),

    # --------------------------------------------------------
    # Concrete stem wall
    #
    # Located beneath the wall/floor junction.
    # --------------------------------------------------------

    RectangleRegion(
        name="concrete_foundation",
        tag=TAG_CONCRETE,
        x=x_ext_ins_1,
        y=-FLOOR_T - FOUNDATION_HEIGHT,
        width=FOUNDATION_T,
        height=FOUNDATION_HEIGHT,
    ),

    # --------------------------------------------------------
    # Exterior foundation insulation
    #
    # Same material tag as wall exterior insulation.
    # --------------------------------------------------------

    RectangleRegion(
        name="foundation_exterior_insulation",
        tag=TAG_EXT_INS,
        x=0.0,
        y=-FLOOR_T - FOUNDATION_HEIGHT,
        width=EXTERIOR_INSULATION_T,
        height=FOUNDATION_HEIGHT + FLOOR_T,
    ),
]

