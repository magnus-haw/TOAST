"""Geometry-independent DOLFINx heat-equation solvers."""

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import ufl
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import VTXWriter

from toast.boundaries import ConvectionBoundary
from toast.materials import Material
from toast.regions import Region


@dataclass(frozen=True)
class SolverConfig:
    initial_temperature: float = 293.15
    dt: float = 60.0
    t_end: float = 24.0 * 3600.0
    volumetric_source: float = 0.0
    degree: int = 1
    output_interval: float | None = 600.0
    output_path: str | Path | None = "wall_temperature.bp"
    petsc_options_prefix: str = "thermal_"


@dataclass
class ThermalResult:
    temperature: fem.Function
    function_space: object
    conductivity: fem.Function
    density: fem.Function
    heat_capacity: fem.Function
    time: float


def build_material_fields(
    domain,
    cell_tags,
    regions: Sequence[Region],
    materials: Mapping[str, Material],
):
    """Create DG0 k, rho and cp fields from semantic region/material mappings."""

    Q0 = fem.functionspace(domain, ("DG", 0))
    k = fem.Function(Q0)
    rho = fem.Function(Q0)
    cp = fem.Function(Q0)
    k.name = "thermal_conductivity"
    rho.name = "density"
    cp.name = "specific_heat"

    k.x.array[:] = np.nan
    rho.x.array[:] = np.nan
    cp.x.array[:] = np.nan

    for region in regions:
        if region.material not in materials:
            raise KeyError(
                f"Region {region.name!r} references unknown material "
                f"{region.material!r}."
            )
        material = materials[region.material]
        cells = cell_tags.find(region.tag)
        k.x.array[cells] = material.k
        rho.x.array[cells] = material.rho
        cp.x.array[cells] = material.cp

    local_bad = (
        np.isnan(k.x.array).any()
        or np.isnan(rho.x.array).any()
        or np.isnan(cp.x.array).any()
    )
    global_bad = domain.comm.allreduce(local_bad, op=MPI.LOR)
    if global_bad:
        raise RuntimeError("At least one mesh cell has no assigned material properties.")

    return k, rho, cp


def _build_steady_forms(domain, facet_tags, V, k, boundaries, volumetric_source=0.0):
    T = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    dx = ufl.Measure("dx", domain=domain)
    ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)

    source = fem.Constant(domain, PETSc.ScalarType(volumetric_source))
    a = k * ufl.dot(ufl.grad(T), ufl.grad(v)) * dx
    L = source * v * dx

    # Keep constants alive until assembly/solve by storing them in a list.
    constants = []
    for bc in boundaries:
        h = fem.Constant(domain, PETSc.ScalarType(bc.h))
        T_inf = fem.Constant(domain, PETSc.ScalarType(bc.temperature))
        constants.extend([h, T_inf])
        a += h * T * v * ds(bc.tag)
        L += h * T_inf * v * ds(bc.tag)

    return a, L, constants


def solve_steady_heat(
    domain,
    cell_tags,
    facet_tags,
    regions: Sequence[Region],
    materials: Mapping[str, Material],
    boundaries: Iterable[ConvectionBoundary],
    *,
    volumetric_source: float = 0.0,
    degree: int = 1,
    petsc_options_prefix: str = "thermal_steady_",
) -> ThermalResult:
    """Solve the steady heterogeneous heat equation with convection boundaries."""

    V = fem.functionspace(domain, ("Lagrange", degree))
    k, rho, cp = build_material_fields(domain, cell_tags, regions, materials)
    T_solution = fem.Function(V)
    T_solution.name = "Temperature"

    a, L, _constants = _build_steady_forms(
        domain,
        facet_tags,
        V,
        k,
        tuple(boundaries),
        volumetric_source,
    )

    problem = LinearProblem(
        a,
        L,
        u=T_solution,
        bcs=[],
        petsc_options_prefix=petsc_options_prefix,
        petsc_options={
            "ksp_type": "cg",
            "pc_type": "gamg",
            "ksp_rtol": 1e-10,
            "ksp_atol": 1e-12,
            "ksp_error_if_not_converged": True,
        },
    )
    problem.solve()
    T_solution.x.scatter_forward()

    return ThermalResult(T_solution, V, k, rho, cp, 0.0)


def solve_transient_heat(
    domain,
    cell_tags,
    facet_tags,
    regions: Sequence[Region],
    materials: Mapping[str, Material],
    boundaries: Iterable[ConvectionBoundary],
    *,
    config: SolverConfig = SolverConfig(),
    progress_every: int | None = 10,
) -> ThermalResult:
    """Backward-Euler transient heat solver for arbitrary tagged 2-D meshes."""

    boundaries = tuple(boundaries)
    V = fem.functionspace(domain, ("Lagrange", config.degree))
    k, rho, cp = build_material_fields(domain, cell_tags, regions, materials)

    T_n = fem.Function(V)
    T_n.name = "Temperature_previous"
    T_solution = fem.Function(V)
    T_solution.name = "Temperature"
    T_n.x.array[:] = config.initial_temperature
    T_solution.x.array[:] = config.initial_temperature

    dt_const = fem.Constant(domain, PETSc.ScalarType(config.dt))
    source = fem.Constant(domain, PETSc.ScalarType(config.volumetric_source))
    T = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    dx = ufl.Measure("dx", domain=domain, subdomain_data=cell_tags)
    ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)

    a = (
        rho * cp / dt_const * T * v * dx
        + k * ufl.dot(ufl.grad(T), ufl.grad(v)) * dx
    )
    L = rho * cp / dt_const * T_n * v * dx + source * v * dx

    constants = [dt_const, source]
    for bc in boundaries:
        h = fem.Constant(domain, PETSc.ScalarType(bc.h))
        T_inf = fem.Constant(domain, PETSc.ScalarType(bc.temperature))
        constants.extend([h, T_inf])
        a += h * T * v * ds(bc.tag)
        L += h * T_inf * v * ds(bc.tag)

    problem = LinearProblem(
        a,
        L,
        u=T_solution,
        bcs=[],
        petsc_options_prefix=config.petsc_options_prefix,
        petsc_options={
            "ksp_type": "cg",
            "pc_type": "gamg",
            "ksp_rtol": 1e-10,
            "ksp_atol": 1e-12,
            "ksp_error_if_not_converged": True,
        },
    )

    writer = None
    if config.output_path is not None:
        writer = VTXWriter(
            domain.comm,
            Path(config.output_path),
            [T_solution],
            engine="BP4",
        )
        writer.write(0.0)

    time = 0.0
    step = 0
    next_output_time = config.output_interval

    try:
        while time < config.t_end - 0.5 * config.dt:
            step += 1
            time += config.dt
            problem.solve()
            T_solution.x.scatter_forward()

            if progress_every and step % progress_every == 0:
                local_min = np.min(T_solution.x.array)
                local_max = np.max(T_solution.x.array)
                Tmin = domain.comm.allreduce(local_min, op=MPI.MIN)
                Tmax = domain.comm.allreduce(local_max, op=MPI.MAX)
                if domain.comm.rank == 0:
                    print(
                        f"step={step:6d} time={time / 3600:8.3f} h "
                        f"Tmin={Tmin - 273.15:8.3f} C "
                        f"Tmax={Tmax - 273.15:8.3f} C"
                    )

            if (
                writer is not None
                and next_output_time is not None
                and time >= next_output_time - 0.5 * config.dt
            ):
                writer.write(time)
                next_output_time += config.output_interval

            T_n.x.array[:] = T_solution.x.array
            T_n.x.scatter_forward()

        if writer is not None:
            writer.write(time)
    finally:
        if writer is not None:
            writer.close()

    # Keep constants referenced through the end of the solve.
    _ = constants
    return ThermalResult(T_solution, V, k, rho, cp, time)
