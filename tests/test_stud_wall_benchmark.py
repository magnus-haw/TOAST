# tests/test_stud_wall_benchmark.py

import numpy as np

from toast.materials import MATERIALS
from toast.regions import Region
from toast.boundaries import ConvectionBoundary

from toast.mesh.stud_wall import (
    StudWallSpec,
    build_stud_wall_mesh,
    EXTERIOR,
    INTERIOR,
)

from toast.solver import (
    solve_steady_heat,
)

from toast.postprocess import (
    boundary_heat_flow,
    boundary_measure,
    effective_u_value,
    energy_balance_error,
    parallel_path_u_value,
    series_u_value,
)

from toast.metrics import convection_boundary_heat_flow


T_EXTERIOR = 273.15
T_INTERIOR = 293.15

H_EXTERIOR = 20.0
H_INTERIOR = 8.0


REGIONS = [
    Region(
        name="insulation",
        material="mineral_wool_ext",
        tag=1,
    ),

    Region(
        name="stud",
        material="wood",
        tag=2,
    ),
]


def test_stud_wall_thermal_bridge():

    spec = StudWallSpec()

    mesh_data = build_stud_wall_mesh(
        spec
    )

    domain = mesh_data.mesh
    cell_tags = mesh_data.cell_tags
    facet_tags = mesh_data.facet_tags

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

    result = solve_steady_heat(
        domain,
        cell_tags,
        facet_tags,
        REGIONS,
        MATERIALS,
        boundaries,
    )

    # -------------------------------------------------------------
    # Integrated heat flow
    # -------------------------------------------------------------

    q_ext = convection_boundary_heat_flow(
        domain,
        result.temperature,
        facet_tags,
        EXTERIOR,
        h=H_EXTERIOR,
        ambient_temperature=T_EXTERIOR,
    )

    q_int = convection_boundary_heat_flow(
        domain,
        result.temperature,
        facet_tags,
        INTERIOR,
        h=H_INTERIOR,
        ambient_temperature=T_INTERIOR,
    )

    # -------------------------------------------------------------
    # Energy conservation
    # -------------------------------------------------------------

    balance = energy_balance_error(
        q_ext,
        q_int,
    )

    # -------------------------------------------------------------
    # Effective U-value
    # -------------------------------------------------------------

    boundary_length = boundary_measure(
        domain,
        facet_tags,
        EXTERIOR,
    )

    delta_T = (
        T_INTERIOR
        - T_EXTERIOR
    )

    U_fem = effective_u_value(
        q_ext,
        boundary_length,
        delta_T,
    )

    # -------------------------------------------------------------
    # Reference homogeneous walls
    # -------------------------------------------------------------

    k_insulation = (
        MATERIALS["mineral_wool_ext"].k
    )

    k_wood = (
        MATERIALS["wood"].k
    )

    U_insulation = series_u_value(
        thickness=spec.thickness,
        conductivity=k_insulation,
        h_exterior=H_EXTERIOR,
        h_interior=H_INTERIOR,
    )

    U_wood = series_u_value(
        thickness=spec.thickness,
        conductivity=k_wood,
        h_exterior=H_EXTERIOR,
        h_interior=H_INTERIOR,
    )

    # -------------------------------------------------------------
    # Parallel path approximation
    # -------------------------------------------------------------

    stud_fraction = (
        spec.stud_width
        / spec.height
    )

    insulation_fraction = (
        1.0
        - stud_fraction
    )

    U_parallel = parallel_path_u_value(
        [
            (
                insulation_fraction,
                k_insulation,
            ),
            (
                stud_fraction,
                k_wood,
            ),
        ],
        thickness=spec.thickness,
        h_exterior=H_EXTERIOR,
        h_interior=H_INTERIOR,
    )

    if domain.comm.rank == 0:

        print()
        print("Stud wall benchmark")
        print("===================")

        print(
            f"Exterior heat flow: "
            f"{q_ext:.6f} W/m"
        )

        print(
            f"Interior heat flow: "
            f"{q_int:.6f} W/m"
        )

        print(
            f"Energy imbalance: "
            f"{balance:.3e}"
        )

        print()

        print(
            f"All-insulation U: "
            f"{U_insulation:.4f} W/(m² K)"
        )

        print(
            f"Parallel-path U:  "
            f"{U_parallel:.4f} W/(m² K)"
        )

        print(
            f"FEM U:            "
            f"{U_fem:.4f} W/(m² K)"
        )

        print(
            f"All-wood U:       "
            f"{U_wood:.4f} W/(m² K)"
        )

    # -------------------------------------------------------------
    # Assertions
    # -------------------------------------------------------------

    # Steady-state energy conservation
    assert balance < 1e-6

    # Stud must increase heat transfer above pure insulation
    assert U_fem > U_insulation

    # But remain below an entirely wood wall
    assert U_fem < U_wood

    # FEM should remain reasonably close to the common engineering
    # parallel-path approximation.
    relative_difference = abs(
        U_fem - U_parallel
    ) / U_parallel

    assert relative_difference < 0.20


