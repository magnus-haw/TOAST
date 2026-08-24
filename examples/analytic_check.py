"""Solve a small layered wall at steady state and compare with the analytic profile."""

from functools import partial

from toast.boundaries import ConvectionBoundary
from toast.materials import Material
from toast.mesh.layered_wall import LayeredWallSpec, build_layered_wall_mesh
from toast.postprocess import analytic_layered_temperature
from toast.regions import Layer
from toast.solver import solve_steady_heat
from toast.visualization import plot_temperature_profile


def main():
    materials = {
        "a": Material(k=0.10, rho=1000.0, cp=1000.0),
        "b": Material(k=0.20, rho=800.0, cp=1200.0),
        "c": Material(k=0.05, rho=40.0, cp=1400.0),
    }
    layers = (
        Layer("a", "a", 1, 0.02, 4),
        Layer("b", "b", 2, 0.03, 4),
        Layer("c", "c", 3, 0.02, 4),
    )
    wall = LayeredWallSpec(layers=layers, height=0.08, vertical_mesh_size=0.02)
    mesh_data = build_layered_wall_mesh(wall, verbose=False)

    T_ext, T_int = 273.15, 293.15
    h_ext, h_int = 20.0, 8.0
    result = solve_steady_heat(
        mesh_data.mesh,
        mesh_data.cell_tags,
        mesh_data.facet_tags,
        wall.layers,
        materials,
        [
            ConvectionBoundary(wall.boundary_tags.exterior, h_ext, T_ext),
            ConvectionBoundary(wall.boundary_tags.interior, h_int, T_int),
        ],
        petsc_options_prefix="analytic_demo_",
    )

    analytic = partial(
        analytic_layered_temperature,
        layers=wall.layers,
        materials=materials,
        h_exterior=h_ext,
        h_interior=h_int,
        T_exterior=T_ext,
        T_interior=T_int,
    )
    plot_temperature_profile(
        mesh_data.mesh,
        result.temperature,
        wall_width=wall.width,
        wall_height=wall.height,
        layers=wall.layers,
        celsius=True,
        analytic=analytic,
        title="Steady layered wall: DOLFINx vs analytic",
    )


if __name__ == "__main__":
    main()
