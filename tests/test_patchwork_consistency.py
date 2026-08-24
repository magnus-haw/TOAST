# tests/test_patchwork_consistency.py

import numpy as np
import ufl

from mpi4py import MPI
from petsc4py import PETSc

from dolfinx import fem, mesh
from dolfinx.fem.petsc import LinearProblem

from toast.mesh.patchwork import (
    RectangleRegion,
    build_rectangular_patchwork_mesh,
)

from toast.verification import (
    manufactured_temperature_ufl,
    manufactured_source_ufl,
)


LX = 1.0
LY = 0.8

K = 0.15

T0 = 293.15
AMPLITUDE = 10.0


REGIONS = [
    RectangleRegion(
        name="bottom_left",
        tag=1,
        x=0.0,
        y=0.0,
        width=LX / 2,
        height=LY / 2,
    ),
    RectangleRegion(
        name="bottom_right",
        tag=2,
        x=LX / 2,
        y=0.0,
        width=LX / 2,
        height=LY / 2,
    ),
    RectangleRegion(
        name="top_left",
        tag=3,
        x=0.0,
        y=LY / 2,
        width=LX / 2,
        height=LY / 2,
    ),
    RectangleRegion(
        name="top_right",
        tag=4,
        x=LX / 2,
        y=LY / 2,
        width=LX / 2,
        height=LY / 2,
    ),
]


def solve_patchwork(mesh_size):

    mesh_data = build_rectangular_patchwork_mesh(
        REGIONS,
        mesh_size=mesh_size,
    )

    domain = mesh_data.mesh
    cell_tags = mesh_data.cell_tags

    V = fem.functionspace(
        domain,
        ("Lagrange", 1),
    )

    Q0 = fem.functionspace(
        domain,
        ("DG", 0),
    )

    # -------------------------------------------------------------
    # Piecewise conductivity field
    #
    # All regions intentionally get the SAME value.
    # -------------------------------------------------------------

    k = fem.Function(Q0)
    k.x.array[:] = np.nan

    for region in REGIONS:

        cells = cell_tags.find(region.tag)

        k.x.array[cells] = K

    assert not np.isnan(k.x.array).any()

    # -------------------------------------------------------------
    # PDE
    # -------------------------------------------------------------

    T = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    Q = manufactured_source_ufl(
        domain,
        lx=LX,
        ly=LY,
        conductivity=K,
        amplitude=AMPLITUDE,
    )

    a = (
        k
        * ufl.dot(
            ufl.grad(T),
            ufl.grad(v),
        )
        * ufl.dx
    )

    L = Q * v * ufl.dx

    # -------------------------------------------------------------
    # Exact boundary value
    #
    # Since the manufactured sine vanishes on all edges,
    # T = T0 everywhere on the boundary.
    # -------------------------------------------------------------

    boundary_facets = mesh.locate_entities_boundary(
        domain,
        domain.topology.dim - 1,
        lambda x: np.full(
            x.shape[1],
            True,
            dtype=bool,
        ),
    )

    boundary_dofs = fem.locate_dofs_topological(
        V,
        domain.topology.dim - 1,
        boundary_facets,
    )

    bc_value = fem.Constant(
        domain,
        PETSc.ScalarType(T0),
    )

    bc = fem.dirichletbc(
        bc_value,
        boundary_dofs,
        V,
    )

    T_h = fem.Function(V)

    problem = LinearProblem(
        a,
        L,
        u=T_h,
        bcs=[bc],
        petsc_options_prefix="patchwork_verify_",
        petsc_options={
            "ksp_type": "preonly",
            "pc_type": "lu",
        },
    )

    problem.solve()

    T_h.x.scatter_forward()

    # -------------------------------------------------------------
    # L2 error
    # -------------------------------------------------------------

    exact = manufactured_temperature_ufl(
        domain,
        lx=LX,
        ly=LY,
        t0=T0,
        amplitude=AMPLITUDE,
    )

    error_form = fem.form(
        (T_h - exact) ** 2
        * ufl.dx
    )

    local_error_sq = fem.assemble_scalar(
        error_form
    )

    global_error_sq = domain.comm.allreduce(
        local_error_sq,
        op=MPI.SUM,
    )

    error = np.sqrt(
        global_error_sq
    )

    return error


def test_patchwork_manufactured_solution_converges():

    error_coarse = solve_patchwork(
        mesh_size=0.08,
    )

    error_fine = solve_patchwork(
        mesh_size=0.04,
    )

    ratio = (
        error_coarse
        / error_fine
    )

    if MPI.COMM_WORLD.rank == 0:

        print()
        print(
            f"Patchwork coarse L2 error: "
            f"{error_coarse:.6e}"
        )

        print(
            f"Patchwork fine L2 error:   "
            f"{error_fine:.6e}"
        )

        print(
            f"Patchwork error ratio:     "
            f"{ratio:.3f}"
        )

    assert error_fine < error_coarse

    assert ratio > 3.0, (
        "Patchwork geometry did not show expected "
        f"P1 L2 convergence; ratio={ratio:.3f}"
    )
