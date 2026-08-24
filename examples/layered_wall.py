"""Run the current layered-wall transient simulation."""

from toast.boundaries import ConvectionBoundary
from toast.materials import MATERIALS
from toast.mesh.layered_wall import DEFAULT_WALL, build_layered_wall_mesh
from toast.solver import SolverConfig, solve_transient_heat
from toast.units import INCH
from toast.visualization import (
    plot_material_regions,
    plot_temperature_field,
    plot_temperature_profile,
)


T_INITIAL = 293.15
T_INTERIOR = 295.15
H_INTERIOR = 8.0
T_EXTERIOR = 273.15
H_EXTERIOR = 20.0

# Turn these off for noninteractive/batch runs.
PLOT_MATERIALS = True
PLOT_TEMPERATURE = True
PLOT_PROFILE = True


def main():
    wall = DEFAULT_WALL
    mesh_data = build_layered_wall_mesh(wall, msh_path="wall.msh")
    domain = mesh_data.mesh

    if domain.comm.rank == 0:
        n_cells = domain.topology.index_map(domain.topology.dim).size_global
        print(f"Mesh loaded\nNumber of cells: {n_cells}")

    boundaries = [
        ConvectionBoundary(
            wall.boundary_tags.exterior,
            H_EXTERIOR,
            T_EXTERIOR,
            "exterior",
        ),
        ConvectionBoundary(
            wall.boundary_tags.interior,
            H_INTERIOR,
            T_INTERIOR,
            "interior",
        ),
    ]

    result = solve_transient_heat(
        domain,
        mesh_data.cell_tags,
        mesh_data.facet_tags,
        wall.layers,
        MATERIALS,
        boundaries,
        config=SolverConfig(
            initial_temperature=T_INITIAL,
            dt=60.0,
            t_end=24.0 * 3600.0,
            output_interval=600.0,
            output_path="wall_temperature.bp",
            petsc_options_prefix="wall_heat_",
        ),
    )

    if domain.comm.rank == 0:
        print("Simulation complete")
        print("Results written to: wall_temperature.bp")

    if PLOT_MATERIALS:
        plot_material_regions(domain, mesh_data.cell_tags)
    if PLOT_TEMPERATURE:
        plot_temperature_field(
            domain,
            result.temperature,
            title=f"Wall temperature after {result.time / 3600:.1f} h",
        )
    if PLOT_PROFILE:
        plot_temperature_profile(
            domain,
            result.temperature,
            wall_width=wall.width,
            wall_height=wall.height,
            layers=wall.layers,
            unit_scale=1.0 / INCH,
            unit_label="in",
            npoints=1000,
            y_fraction=0.5,
            title=f"Mid-height temperature profile at {result.time / 3600:.1f} h",
        )


if __name__ == "__main__":
    main()
