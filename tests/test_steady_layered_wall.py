import numpy as np

from toast.boundaries import ConvectionBoundary
from toast.materials import Material
from toast.mesh.layered_wall import LayeredWallSpec, build_layered_wall_mesh
from toast.postprocess import analytic_layered_temperature, sample_horizontal_profile
from toast.regions import Layer
from toast.solver import solve_steady_heat


def test_steady_layered_wall_matches_analytic_solution():
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
    wall = LayeredWallSpec(
        layers=layers,
        height=0.08,
        vertical_mesh_size=0.02,
    )
    mesh_data = build_layered_wall_mesh(wall, verbose=False)

    T_ext = 273.15
    T_int = 293.15
    h_ext = 20.0
    h_int = 8.0
    boundaries = [
        ConvectionBoundary(wall.boundary_tags.exterior, h_ext, T_ext),
        ConvectionBoundary(wall.boundary_tags.interior, h_int, T_int),
    ]

    result = solve_steady_heat(
        mesh_data.mesh,
        mesh_data.cell_tags,
        mesh_data.facet_tags,
        wall.layers,
        materials,
        boundaries,
        petsc_options_prefix="test_steady_",
    )

    x, T_numeric, valid = sample_horizontal_profile(
        mesh_data.mesh,
        result.temperature,
        width=wall.width,
        height=wall.height,
        npoints=101,
    )
    T_exact = analytic_layered_temperature(
        x,
        wall.layers,
        materials,
        h_ext,
        h_int,
        T_ext,
        T_int,
    )

    error = np.max(np.abs(T_numeric[valid] - T_exact[valid]))
    assert error < 0.10
