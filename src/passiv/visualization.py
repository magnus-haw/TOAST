"""Visualization helpers adapted from the working visualization.py."""

import numpy as np
import matplotlib.pyplot as plt

from toast.postprocess import sample_horizontal_profile


def plot_temperature_field(
    domain,
    temperature,
    *,
    title="Temperature field",
    show_edges=False,
    celsius=True,
):
    """Interactive 2-D PyVista temperature field for CG1 temperature meshes."""

    import pyvista
    from dolfinx import plot

    topology, cell_types, coordinates = plot.vtk_mesh(
        domain, domain.topology.dim
    )
    grid = pyvista.UnstructuredGrid(topology, cell_types, coordinates)
    values = temperature.x.array.real.copy()

    if celsius:
        values -= 273.15
        field_name = "Temperature [C]"
        scalar_title = "Temperature [°C]"
    else:
        field_name = "Temperature [K]"
        scalar_title = "Temperature [K]"

    if len(values) != grid.n_points:
        raise RuntimeError(
            f"Temperature DOFs ({len(values)}) do not match "
            f"mesh VTK points ({grid.n_points})."
        )

    grid.point_data[field_name] = values
    plotter = pyvista.Plotter()
    plotter.add_mesh(
        grid,
        scalars=field_name,
        show_edges=show_edges,
        cmap="coolwarm",
        scalar_bar_args={"title": scalar_title},
    )
    plotter.add_title(title)
    plotter.view_xy()
    plotter.show_axes()
    plotter.show()


def plot_material_regions(
    domain,
    cell_tags,
    *,
    title="Material regions",
    show_edges=True,
):
    """Interactive PyVista display of integer cell/material tags."""

    import pyvista
    from dolfinx import plot

    tdim = domain.topology.dim
    topology, cell_types, coordinates = plot.vtk_mesh(domain, tdim)
    grid = pyvista.UnstructuredGrid(topology, cell_types, coordinates)
    tags = np.zeros(grid.n_cells, dtype=np.int32)
    valid = cell_tags.indices < grid.n_cells
    tags[cell_tags.indices[valid]] = cell_tags.values[valid]
    grid.cell_data["Material"] = tags

    plotter = pyvista.Plotter()
    plotter.add_mesh(
        grid,
        scalars="Material",
        categories=True,
        show_edges=show_edges,
    )
    plotter.add_title(title)
    plotter.view_xy()
    plotter.show_axes()
    plotter.show()


def plot_temperature_profile(
    domain,
    temperature,
    *,
    wall_width,
    wall_height,
    layers=None,
    unit_scale=1.0,
    unit_label="m",
    npoints=500,
    y_fraction=0.5,
    celsius=True,
    title="Wall temperature profile",
    analytic=None,
):
    """Matplotlib through-wall profile with optional analytic comparison."""

    if not 0.0 <= y_fraction <= 1.0:
        raise ValueError("y_fraction must lie between 0 and 1.")

    x, values, valid = sample_horizontal_profile(
        domain,
        temperature,
        width=wall_width,
        height=wall_height,
        y_fraction=y_fraction,
        npoints=npoints,
    )
    if not np.any(valid):
        raise RuntimeError("Could not locate any profile points inside the mesh.")

    if celsius:
        values = values - 273.15
        ylabel = "Temperature [°C]"
    else:
        ylabel = "Temperature [K]"

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x[valid] * unit_scale, values[valid], linewidth=2, label="DOLFINx")

    if analytic is not None:
        analytic_values = np.asarray(analytic(x), dtype=float)
        if celsius:
            analytic_values -= 273.15
        ax.plot(
            x * unit_scale,
            analytic_values,
            linestyle="--",
            linewidth=1.5,
            label="Analytic",
        )
        ax.legend()

    if layers is not None:
        xpos = 0.0
        for layer in list(layers)[:-1]:
            xpos += layer.thickness
            ax.axvline(xpos * unit_scale, linestyle="--", linewidth=1, alpha=0.3)

    ax.set_xlabel(f"Distance from exterior [{unit_label}]")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    plt.show()
    return x, values


def plot_mesh(domain, *, title="Finite-element mesh", show_edges=True):
    """Interactive PyVista mesh-only display."""

    import pyvista
    from dolfinx import plot

    topology, cell_types, coordinates = plot.vtk_mesh(
        domain, domain.topology.dim
    )
    grid = pyvista.UnstructuredGrid(topology, cell_types, coordinates)
    plotter = pyvista.Plotter()
    plotter.add_mesh(grid, show_edges=show_edges, color="white")
    plotter.add_title(title)
    plotter.view_xy()
    plotter.show_axes()
    plotter.show()


def plot_multiple_temperature_profiles(
    domain,
    temperature,
    *,
    wall_width,
    y_positions,
    unit_scale=1.0,
    unit_label="m",
    npoints=400,
):
    import matplotlib.pyplot as plt

    from .visualization import (
        evaluate_function_at_points,
    )

    x = np.linspace(
        1e-9,
        wall_width - 1e-9,
        npoints,
    )

    fig, ax = plt.subplots()

    for y in y_positions:

        points = np.column_stack([
            x,
            np.full_like(x, y),
        ])

        values, valid = (
            evaluate_function_at_points(
                domain,
                temperature,
                points,
            )
        )

        ax.plot(
            x[valid] * unit_scale,
            values[valid] - 273.15,
            label=f"y={y * unit_scale:.2f} {unit_label}",
        )

    ax.set_xlabel(
        f"Through-wall distance [{unit_label}]"
    )

    ax.set_ylabel(
        "Temperature [°C]"
    )

    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    plt.show()



