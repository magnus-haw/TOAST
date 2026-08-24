import numpy as np

from toast.boundaries import ConvectionBoundary
from toast.materials import Material
from toast.mesh.layered_wall import LayeredWallSpec, build_layered_wall_mesh
from toast.regions import Layer
from toast.solver import SolverConfig, solve_transient_heat


def test_uniform_equilibrium_is_preserved():
    T0 = 293.15
    layers = (Layer("single", "solid", 1, 0.05, 3),)
    materials = {"solid": Material(k=0.5, rho=1000.0, cp=1000.0)}
    wall = LayeredWallSpec(
        layers=layers,
        height=0.05,
        vertical_mesh_size=0.01,
    )
    mesh_data = build_layered_wall_mesh(wall, verbose=False)
    boundaries = [
        ConvectionBoundary(wall.boundary_tags.exterior, 10.0, T0),
        ConvectionBoundary(wall.boundary_tags.interior, 10.0, T0),
    ]

    result = solve_transient_heat(
        mesh_data.mesh,
        mesh_data.cell_tags,
        mesh_data.facet_tags,
        wall.layers,
        materials,
        boundaries,
        config=SolverConfig(
            initial_temperature=T0,
            dt=10.0,
            t_end=100.0,
            output_interval=None,
            output_path=None,
            petsc_options_prefix="test_uniform_",
        ),
        progress_every=None,
    )

    assert np.max(np.abs(result.temperature.x.array - T0)) < 1e-6
