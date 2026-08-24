"""
Thermal-performance metrics and integrated heat-flow utilities.

The functions in this module are geometry-independent and operate on
DOLFINx solutions, tagged boundaries, or scalar reference values.

Sign convention
---------------
boundary_heat_flow() reports conductive heat flow OUTWARD from the
computational domain.

For a steady problem with no volumetric source:

    Q_exterior + Q_interior + ... = 0

For a 2-D cross-sectional simulation, integrated heat flow has units
W/m, where the implicit out-of-plane depth is one meter.
"""

from dataclasses import dataclass

import numpy as np
import ufl

from mpi4py import MPI
from dolfinx import fem


# =====================================================================
# Result containers
# =====================================================================

@dataclass(frozen=True)
class AssemblyMetrics:
    """
    Summary metrics for a steady thermal assembly.

    Attributes
    ----------
    exterior_heat_flow
        Conductive heat flow outward through exterior boundary [W/m].

    interior_heat_flow
        Conductive heat flow outward through interior boundary [W/m].

    energy_balance_error
        Dimensionless relative heat-flow imbalance.

    boundary_length
        Reference boundary length [m].

    delta_temperature
        Interior minus exterior ambient temperature [K].

    u_value
        Effective assembly U-value [W/(m^2 K)].
    """

    exterior_heat_flow: float

    interior_heat_flow: float

    energy_balance_error: float

    boundary_length: float

    delta_temperature: float

    u_value: float


# =====================================================================
# Integrated FEM quantities
# =====================================================================

def boundary_heat_flow(
    domain,
    temperature,
    conductivity,
    facet_tags,
    boundary_tag,
):
    """
    Integrate conductive heat flow outward through a boundary.

    Parameters
    ----------
    domain
        DOLFINx mesh.

    temperature
        DOLFINx temperature Function.

    conductivity
        DOLFINx scalar conductivity Function or compatible UFL
        coefficient.

    facet_tags
        DOLFINx facet MeshTags.

    boundary_tag
        Integer physical tag identifying the boundary.

    Returns
    -------
    float
        Integrated outward conductive heat flow.

        For a 2-D model:
            [W/m]

        assuming one meter of out-of-plane depth.

    Notes
    -----
    Fourier heat flux is

        q = -k grad(T)

    and the outward normal component is

        q_n = -k grad(T) . n.
    """

    n = ufl.FacetNormal(
        domain
    )

    ds = ufl.Measure(
        "ds",
        domain=domain,
        subdomain_data=facet_tags,
    )

    outward_flux = (
        -conductivity
        * ufl.dot(
            ufl.grad(
                temperature
            ),
            n,
        )
    )

    form = fem.form(
        outward_flux
        * ds(boundary_tag)
    )

    local_value = (
        fem.assemble_scalar(
            form
        )
    )

    return float(
        domain.comm.allreduce(
            local_value,
            op=MPI.SUM,
        )
    )


def boundary_measure(
    domain,
    facet_tags,
    boundary_tag,
):
    """
    Integrate the measure of a tagged boundary.

    In a 2-D mesh this returns boundary length [m].
    """

    ds = ufl.Measure(
        "ds",
        domain=domain,
        subdomain_data=facet_tags,
    )

    one = fem.Constant(
        domain,
        1.0,
    )

    form = fem.form(
        one
        * ds(boundary_tag)
    )

    local_value = (
        fem.assemble_scalar(
            form
        )
    )

    return float(
        domain.comm.allreduce(
            local_value,
            op=MPI.SUM,
        )
    )


# =====================================================================
# Energy conservation
# =====================================================================

def energy_balance_error(
    *heat_flows,
):
    """
    Compute relative steady-state heat-flow imbalance.

    Parameters
    ----------
    heat_flows
        Signed outward heat flows through all relevant boundaries.

    Returns
    -------
    float
        Dimensionless relative imbalance:

            abs(sum(Q_i))
            ----------------
            sum(abs(Q_i))

    A perfectly conservative steady solution with zero volumetric
    source gives zero.
    """

    values = np.asarray(
        heat_flows,
        dtype=float,
    )

    if values.size == 0:
        raise ValueError(
            "At least one heat flow is required."
        )

    numerator = abs(
        np.sum(values)
    )

    denominator = np.sum(
        np.abs(values)
    )

    if denominator == 0.0:
        return 0.0

    return float(
        numerator
        / denominator
    )


# =====================================================================
# Effective U-value
# =====================================================================

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
        Integrated heat flow [W/m] from a 2-D model.

    boundary_length
        Reference wall boundary length [m].

    delta_temperature
        Indoor-outdoor temperature difference [K].

    Returns
    -------
    float
        Effective U-value [W/(m^2 K)].
    """

    if boundary_length <= 0.0:
        raise ValueError(
            "boundary_length must be positive."
        )

    if np.isclose(
        delta_temperature,
        0.0,
    ):
        raise ValueError(
            "delta_temperature must be nonzero."
        )

    return float(
        abs(heat_flow)
        / (
            boundary_length
            * abs(delta_temperature)
        )
    )


# =====================================================================
# 1-D analytical references
# =====================================================================

def series_thermal_resistance(
    *,
    thickness,
    conductivity,
    h_exterior=None,
    h_interior=None,
):
    """
    Thermal resistance of one homogeneous layer with optional films.

    Returns
    -------
    float
        Resistance [m^2 K/W].
    """

    if thickness < 0.0:
        raise ValueError(
            "thickness cannot be negative."
        )

    if conductivity <= 0.0:
        raise ValueError(
            "conductivity must be positive."
        )

    resistance = (
        thickness
        / conductivity
    )

    if h_exterior is not None:

        if h_exterior <= 0.0:
            raise ValueError(
                "h_exterior must be positive."
            )

        resistance += (
            1.0
            / h_exterior
        )

    if h_interior is not None:

        if h_interior <= 0.0:
            raise ValueError(
                "h_interior must be positive."
            )

        resistance += (
            1.0
            / h_interior
        )

    return float(
        resistance
    )


def series_u_value(
    *,
    thickness,
    conductivity,
    h_exterior=None,
    h_interior=None,
):
    """
    U-value for a homogeneous one-dimensional wall path.

    Returns
    -------
    float
        U-value [W/(m^2 K)].
    """

    resistance = (
        series_thermal_resistance(
            thickness=thickness,
            conductivity=conductivity,
            h_exterior=h_exterior,
            h_interior=h_interior,
        )
    )

    return float(
        1.0
        / resistance
    )


def multilayer_thermal_resistance(
    layers,
    *,
    h_exterior=None,
    h_interior=None,
):
    """
    Series resistance of multiple planar layers.

    Parameters
    ----------
    layers
        Iterable of ``(thickness, conductivity)``.

    Returns
    -------
    float
        Total thermal resistance [m^2 K/W].
    """

    layers = list(
        layers
    )

    if not layers:
        raise ValueError(
            "At least one layer is required."
        )

    resistance = 0.0

    if h_exterior is not None:

        if h_exterior <= 0.0:
            raise ValueError(
                "h_exterior must be positive."
            )

        resistance += (
            1.0
            / h_exterior
        )

    for (
        thickness,
        conductivity,
    ) in layers:

        if thickness < 0.0:
            raise ValueError(
                "Layer thickness cannot "
                "be negative."
            )

        if conductivity <= 0.0:
            raise ValueError(
                "Layer conductivity must "
                "be positive."
            )

        resistance += (
            thickness
            / conductivity
        )

    if h_interior is not None:

        if h_interior <= 0.0:
            raise ValueError(
                "h_interior must be positive."
            )

        resistance += (
            1.0
            / h_interior
        )

    return float(
        resistance
    )


def multilayer_u_value(
    layers,
    *,
    h_exterior=None,
    h_interior=None,
):
    """
    U-value of multiple one-dimensional layers in series.
    """

    return float(
        1.0
        / multilayer_thermal_resistance(
            layers,
            h_exterior=h_exterior,
            h_interior=h_interior,
        )
    )


# =====================================================================
# Parallel path approximation
# =====================================================================

def parallel_path_u_value(
    paths,
    *,
    thickness,
    h_exterior=None,
    h_interior=None,
):
    """
    Area-weighted parallel-path U-value.

    Parameters
    ----------
    paths
        Iterable of:

            (area_fraction, conductivity)

        For example:

            [
                (0.90, k_insulation),
                (0.10, k_wood),
            ]

    thickness
        Common path thickness [m].

    h_exterior, h_interior
        Optional surface-film coefficients [W/(m^2 K)].

    Returns
    -------
    float
        Parallel-path U-value [W/(m^2 K)].
    """

    paths = list(
        paths
    )

    if not paths:
        raise ValueError(
            "At least one parallel path is required."
        )

    fraction_sum = sum(
        fraction
        for fraction, _
        in paths
    )

    if not np.isclose(
        fraction_sum,
        1.0,
        rtol=0.0,
        atol=1e-10,
    ):
        raise ValueError(
            "Parallel-path fractions must "
            "sum to 1. "
            f"Received {fraction_sum:.12g}."
        )

    U = 0.0

    for (
        fraction,
        conductivity,
    ) in paths:

        if fraction < 0.0:
            raise ValueError(
                "Path fractions cannot "
                "be negative."
            )

        U_path = series_u_value(
            thickness=thickness,
            conductivity=conductivity,
            h_exterior=h_exterior,
            h_interior=h_interior,
        )

        U += (
            fraction
            * U_path
        )

    return float(U)


# =====================================================================
# 2-D junction / thermal bridge metrics
# =====================================================================

def thermal_coupling_coefficient(
    heat_flow,
    delta_temperature,
):
    """
    Thermal coupling coefficient from a 2-D calculation.

    Returns
    -------
    float
        [W/(m K)] for a 2-D model with unit out-of-plane depth.
    """

    if np.isclose(
        delta_temperature,
        0.0,
    ):
        raise ValueError(
            "delta_temperature must be nonzero."
        )

    return float(
        abs(heat_flow)
        / abs(delta_temperature)
    )


def linear_thermal_transmittance(
    heat_flow,
    delta_temperature,
    reference_u_lengths,
):
    """
    Compute linear thermal transmittance psi.

    Parameters
    ----------
    heat_flow
        Total 2-D heat flow [W/m].

    delta_temperature
        Indoor-outdoor temperature difference [K].

    reference_u_lengths
        Iterable of:

            (U_i, L_i)

        where U_i is a reference planar U-value [W/(m^2 K)]
        and L_i is its associated model length [m].

    Returns
    -------
    float
        Linear thermal transmittance psi [W/(m K)].

    Notes
    -----
    psi = L_2D - sum(U_i L_i)

    where

        L_2D = Q / DeltaT.
    """

    reference_u_lengths = list(
        reference_u_lengths
    )

    if not reference_u_lengths:
        raise ValueError(
            "At least one reference "
            "(U, length) term is required."
        )

    coupling = (
        thermal_coupling_coefficient(
            heat_flow,
            delta_temperature,
        )
    )

    one_d_reference = 0.0

    for (
        U,
        length,
    ) in reference_u_lengths:

        if U < 0.0:
            raise ValueError(
                "Reference U-values cannot "
                "be negative."
            )

        if length < 0.0:
            raise ValueError(
                "Reference lengths cannot "
                "be negative."
            )

        one_d_reference += (
            U
            * length
        )

    return float(
        coupling
        - one_d_reference
    )


# =====================================================================
# Convenience assembly calculation
# =====================================================================

def compute_assembly_metrics(
    *,
    domain,
    temperature,
    conductivity,
    facet_tags,
    exterior_tag,
    interior_tag,
    exterior_temperature,
    interior_temperature,
    h_exterior,
    h_interior,
):
    """
    Compute common steady-state assembly metrics.

    Parameters
    ----------
    domain
        DOLFINx mesh.

    temperature
        Solved DOLFINx temperature Function.

    conductivity
        DOLFINx conductivity field.

    facet_tags
        Boundary MeshTags.

    exterior_tag, interior_tag
        Boundary tags.

    exterior_temperature, interior_temperature
        Ambient temperatures [K].

    Returns
    -------
    AssemblyMetrics
    """

    q_exterior = convection_boundary_heat_flow(
        domain,
        temperature,
        facet_tags,
        exterior_tag,
        h=h_exterior,
        ambient_temperature=exterior_temperature,
    )

    q_interior = convection_boundary_heat_flow(
        domain,
        temperature,
        facet_tags,
        interior_tag,
        h=h_interior,
        ambient_temperature=interior_temperature,
    )

    length_exterior = boundary_measure(
        domain,
        facet_tags,
        exterior_tag,
    )

    length_interior = boundary_measure(
        domain,
        facet_tags,
        interior_tag,
    )

    # For an ordinary wall cross section these should be identical.
    if not np.isclose(
        length_exterior,
        length_interior,
        rtol=1e-8,
        atol=1e-12,
    ):
        raise ValueError(
            "Exterior and interior reference "
            "boundary lengths differ: "
            f"{length_exterior} vs "
            f"{length_interior}. "
            "For a junction geometry, compute "
            "the reference metric explicitly."
        )

    delta_temperature = (
        interior_temperature
        - exterior_temperature
    )

    balance = energy_balance_error(
        q_exterior,
        q_interior,
    )

    U = effective_u_value(
        q_exterior,
        length_exterior,
        delta_temperature,
    )

    return AssemblyMetrics(
        exterior_heat_flow=q_exterior,
        interior_heat_flow=q_interior,
        energy_balance_error=balance,
        boundary_length=length_exterior,
        delta_temperature=float(
            delta_temperature
        ),
        u_value=U,
    )

def convection_boundary_heat_flow(
    domain,
    temperature,
    facet_tags,
    boundary_tag,
    *,
    h,
    ambient_temperature,
):
    """
    Integrate outward heat flow through a convection boundary.

    Boundary condition
    ------------------
        -k grad(T) . n = h (T - T_inf)

    Therefore the outward heat flux is evaluated directly as

        q_out = h (T - T_inf).

    Parameters
    ----------
    domain
        DOLFINx mesh.

    temperature
        Solved temperature Function [K].

    facet_tags
        Boundary MeshTags.

    boundary_tag
        Integer facet tag.

    h
        Convection coefficient [W/(m^2 K)].

    ambient_temperature
        Ambient temperature T_inf [K].

    Returns
    -------
    float
        Integrated outward heat flow [W/m] for a 2-D model.
    """

    ds = ufl.Measure(
        "ds",
        domain=domain,
        subdomain_data=facet_tags,
    )

    q_out = (
        h
        * (
            temperature
            - ambient_temperature
        )
    )

    form = fem.form(
        q_out
        * ds(boundary_tag)
    )

    local_value = fem.assemble_scalar(
        form
    )

    return float(
        domain.comm.allreduce(
            local_value,
            op=MPI.SUM,
        )
    )


