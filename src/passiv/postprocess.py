# src/toast/postprocess.py
"""Post-processing, sampling and analytic layered-wall verification tools."""

import numpy as np
from dolfinx import geometry
import ufl

from mpi4py import MPI
from dolfinx import fem


def kelvin_to_celsius(value):
    return np.asarray(value) - 273.15


def layered_wall_resistance(layers, materials, h_exterior, h_interior):
    """Total area-normalized resistance [m^2 K/W] for 1-D layers in series."""

    return (
        1.0 / h_exterior
        + sum(layer.thickness / materials[layer.material].k for layer in layers)
        + 1.0 / h_interior
    )


def layered_wall_heat_flux(
    layers,
    materials,
    h_exterior,
    h_interior,
    T_exterior,
    T_interior,
):
    """Positive heat flux from interior toward exterior [W/m^2]."""

    resistance = layered_wall_resistance(
        layers, materials, h_exterior, h_interior
    )
    return (T_interior - T_exterior) / resistance


def analytic_layered_temperature(
    x,
    layers,
    materials,
    h_exterior,
    h_interior,
    T_exterior,
    T_interior,
):
    """Piecewise-linear steady temperature inside a layered wall.

    Coordinate x starts at the exterior material surface and increases toward
    the interior. Convective film drops are included in the heat flux and in
    the exterior-surface starting temperature.
    """

    x = np.asarray(x, dtype=float)
    width = sum(layer.thickness for layer in layers)
    if np.any(x < -1e-12) or np.any(x > width + 1e-12):
        raise ValueError("x must lie within the wall thickness.")

    q = layered_wall_heat_flux(
        layers,
        materials,
        h_exterior,
        h_interior,
        T_exterior,
        T_interior,
    )
    T_surface_ext = T_exterior + q / h_exterior

    values = np.empty_like(x)
    x0 = 0.0
    T0 = T_surface_ext
    for i, layer in enumerate(layers):
        x1 = x0 + layer.thickness
        if i == len(layers) - 1:
            mask = (x >= x0 - 1e-14) & (x <= x1 + 1e-14)
        else:
            mask = (x >= x0 - 1e-14) & (x < x1)
        k = materials[layer.material].k
        values[mask] = T0 + q * (x[mask] - x0) / k
        T0 += q * layer.thickness / k
        x0 = x1

    return values


def evaluate_function_at_points(domain, function, points):
    """Evaluate a scalar DOLFINx function at N physical points.

    Intended primarily for serial post-processing. Points may be (N, 2) or
    (N, 3); 2-D points are padded to 3-D coordinates for DOLFINx geometry calls.
    """

    points = np.asarray(points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] not in (2, 3):
        raise ValueError("points must have shape (N, 2) or (N, 3)")
    if points.shape[1] == 2:
        padded = np.zeros((points.shape[0], 3), dtype=np.float64)
        padded[:, :2] = points
        points = padded

    tree = geometry.bb_tree(domain, domain.topology.dim)
    candidates = geometry.compute_collisions_points(tree, points)
    colliding_cells = geometry.compute_colliding_cells(
        domain, candidates, points
    )

    values = np.full(points.shape[0], np.nan, dtype=np.float64)
    valid = np.zeros(points.shape[0], dtype=bool)
    for i, point in enumerate(points):
        cells = colliding_cells.links(i)
        if len(cells) == 0:
            continue
        value = function.eval(point, cells[0])
        values[i] = np.asarray(value).flat[0].real
        valid[i] = True
    return values, valid


def sample_horizontal_profile(
    domain,
    function,
    *,
    width,
    height,
    y_fraction=0.5,
    npoints=500,
):
    """Sample a scalar field along x at a fixed fraction of the geometry height."""

    eps = max(1e-10, width * 1e-10)
    x = np.linspace(eps, width - eps, npoints)
    points = np.column_stack([x, np.full_like(x, y_fraction * height)])
    values, valid = evaluate_function_at_points(domain, function, points)
    return x, values, valid




def boundary_heat_flow(
    domain,
    temperature,
    conductivity,
    facet_tags,
    boundary_tag,
):
    """
    Integrate conductive heat flow OUTWARD through a tagged boundary.

    For a 2-D model, the result has units W/m, assuming unit depth
    perpendicular to the modeled plane.

    Sign convention
    ---------------
    Positive:
        net heat leaving the computational domain.

    Negative:
        net heat entering the computational domain.
    """

    n = ufl.FacetNormal(domain)

    ds = ufl.Measure(
        "ds",
        domain=domain,
        subdomain_data=facet_tags,
    )

    outward_flux = (
        -conductivity
        * ufl.dot(
            ufl.grad(temperature),
            n,
        )
    )

    form = fem.form(
        outward_flux * ds(boundary_tag)
    )

    local_value = fem.assemble_scalar(form)

    return domain.comm.allreduce(
        local_value,
        op=MPI.SUM,
    )


def boundary_measure(
    domain,
    facet_tags,
    boundary_tag,
):
    """
    Return the length of a tagged boundary in a 2-D model.

    For a unit-depth interpretation this is numerically equal to the
    boundary area in m^2 per meter of out-of-plane depth.
    """

    ds = ufl.Measure(
        "ds",
        domain=domain,
        subdomain_data=facet_tags,
    )

    form = fem.form(
        fem.Constant(domain, 1.0)
        * ds(boundary_tag)
    )

    local_value = fem.assemble_scalar(form)

    return domain.comm.allreduce(
        local_value,
        op=MPI.SUM,
    )


def effective_u_value(
    heat_flow,
    boundary_length,
    delta_temperature,
):
    """
    Compute effective assembly U-value.

    Parameters
    ----------
    heat_flow
        Total heat flow per unit out-of-plane depth [W/m].

    boundary_length
        Interior/exterior boundary length [m].

    delta_temperature
        Indoor-outdoor air-temperature difference [K].

    Returns
    -------
    float
        U-value [W/(m^2 K)].
    """

    if boundary_length <= 0:
        raise ValueError(
            "boundary_length must be positive."
        )

    if np.isclose(delta_temperature, 0.0):
        raise ValueError(
            "delta_temperature must be nonzero."
        )

    return (
        abs(heat_flow)
        / (
            boundary_length
            * abs(delta_temperature)
        )
    )


def energy_balance_error(
    *heat_flows,
):
    """
    Return relative steady-state boundary energy imbalance.

    For zero volumetric heating at steady state,

        sum(Q_out) = 0.

    The returned value is dimensionless.
    """

    heat_flows = np.asarray(
        heat_flows,
        dtype=float,
    )

    total = abs(
        np.sum(heat_flows)
    )

    scale = np.sum(
        np.abs(heat_flows)
    )

    if scale == 0:
        return 0.0

    return total / scale


def series_u_value(
    *,
    thickness,
    conductivity,
    h_exterior,
    h_interior,
):
    """
    One-dimensional U-value for a homogeneous wall.
    """

    resistance = (
        1.0 / h_exterior
        + thickness / conductivity
        + 1.0 / h_interior
    )

    return 1.0 / resistance


def parallel_path_u_value(
    paths,
    *,
    thickness,
    h_exterior,
    h_interior,
):
    """
    Calculate a simple area-weighted parallel-path U-value.

    Parameters
    ----------
    paths
        Iterable of (fraction, conductivity).

        Fractions should sum to 1.

    thickness
        Wall thickness [m].
    """

    paths = list(paths)

    fraction_sum = sum(
        fraction
        for fraction, _ in paths
    )

    if not np.isclose(
        fraction_sum,
        1.0,
        atol=1e-12,
    ):
        raise ValueError(
            f"Path fractions sum to {fraction_sum}, not 1."
        )

    U = 0.0

    for fraction, conductivity in paths:

        U_path = series_u_value(
            thickness=thickness,
            conductivity=conductivity,
            h_exterior=h_exterior,
            h_interior=h_interior,
        )

        U += fraction * U_path

    return U