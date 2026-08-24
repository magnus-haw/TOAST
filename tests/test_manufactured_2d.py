"""
Verify the 2-D steady heat solver against a manufactured analytic solution.

PDE
---
    -div(k grad(T)) = Q

Exact solution
--------------
    T(x, y) = T0
              + A sin(pi*x/Lx)
                  sin(pi*y/Ly)

with T = T0 on the entire boundary.

For linear P1 finite elements, the L2 error should converge approximately
as O(h^2).
"""

import numpy as np
import ufl

from mpi4py import MPI
from petsc4py import PETSc

from dolfinx import fem, mesh
from dolfinx.fem.petsc import LinearProblem


# ----------------------------------------------------------------------
# Problem definition
# ----------------------------------------------------------------------

LX = 1.0
LY = 0.7

K = 0.15

T0 = 293.15
AMPLITUDE = 10.0


# ----------------------------------------------------------------------
# Exact solution
# ----------------------------------------------------------------------

def exact_temperature_numpy(x):
    """
    NumPy version used for interpolation/evaluation.
    """

    return (
        T0
        + AMPLITUDE
        * np.sin(np.pi * x[0] / LX)
        * np.sin(np.pi * x[1] / LY)
    )


# ----------------------------------------------------------------------
# Solve one mesh resolution
# ----------------------------------------------------------------------

def solve_manufactured_problem(nx, ny):

    domain = mesh.create_rectangle(
        MPI.COMM_WORLD,
        [
            np.array([0.0, 0.0]),
            np.array([LX, LY]),
        ],
        [nx, ny],
        cell_type=mesh.CellType.triangle,
    )

    V = fem.functionspace(
        domain,
        ("Lagrange", 1),
    )

    # --------------------------------------------------------------
    # Trial/test functions
    # --------------------------------------------------------------

    T = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    x = ufl.SpatialCoordinate(domain)

    # --------------------------------------------------------------
    # Manufactured source
    # --------------------------------------------------------------

    spatial_shape = (
        ufl.sin(np.pi * x[0] / LX)
        * ufl.sin(np.pi * x[1] / LY)
    )

    Q = (
        K
        * AMPLITUDE
        * np.pi**2
        * (
            1.0 / LX**2
            + 1.0 / LY**2
        )
        * spatial_shape
    )

    # --------------------------------------------------------------
    # Weak form
    #
    # -div(k grad T) = Q
    #
    # becomes
    #
    # int k grad(T).grad(v) dx
    # =
    # int Q v dx
    # --------------------------------------------------------------

    a = (
        K
        * ufl.dot(
            ufl.grad(T),
            ufl.grad(v),
        )
        * ufl.dx
    )

    L = Q * v * ufl.dx

    # --------------------------------------------------------------
    # Dirichlet BC: T = T0 on entire boundary
    # --------------------------------------------------------------

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

    T_boundary = fem.Constant(
        domain,
        PETSc.ScalarType(T0),
    )

    bc = fem.dirichletbc(
        T_boundary,
        boundary_dofs,
        V,
    )

    # --------------------------------------------------------------
    # Solve
    # --------------------------------------------------------------

    T_h = fem.Function(V)
    T_h.name = "Temperature"

    problem = LinearProblem(
        a,
        L,
        bcs=[bc],
        u=T_h,
        petsc_options_prefix="manufactured_heat_",
        petsc_options={
            "ksp_type": "preonly",
            "pc_type": "lu",
            "ksp_error_if_not_converged": True,
        },
    )

    problem.solve()

    T_h.x.scatter_forward()

    # --------------------------------------------------------------
    # Exact solution represented symbolically for error integration
    # --------------------------------------------------------------

    T_exact_ufl = (
        T0
        + AMPLITUDE
        * ufl.sin(np.pi * x[0] / LX)
        * ufl.sin(np.pi * x[1] / LY)
    )

    # --------------------------------------------------------------
    # L2 error
    # --------------------------------------------------------------

    error_form = fem.form(
        (T_h - T_exact_ufl) ** 2
        * ufl.dx
    )

    local_error_sq = fem.assemble_scalar(
        error_form
    )

    global_error_sq = domain.comm.allreduce(
        local_error_sq,
        op=MPI.SUM,
    )

    l2_error = np.sqrt(
        global_error_sq
    )

    return domain, V, T_h, l2_error


# ----------------------------------------------------------------------
# Verification test
# ----------------------------------------------------------------------

def test_manufactured_solution_converges():

    # Coarse mesh
    _, _, _, error_coarse = (
        solve_manufactured_problem(
            nx=12,
            ny=10,
        )
    )

    # Refine each direction by factor 2
    _, _, _, error_fine = (
        solve_manufactured_problem(
            nx=24,
            ny=20,
        )
    )

    ratio = error_coarse / error_fine

    if MPI.COMM_WORLD.rank == 0:
        print()
        print(
            f"Coarse L2 error : {error_coarse:.6e}"
        )
        print(
            f"Fine L2 error   : {error_fine:.6e}"
        )
        print(
            f"Error ratio     : {ratio:.3f}"
        )

    # P1 elements should approach second-order L2 convergence.
    #
    # Perfect asymptotic behavior would give:
    #
    # error(h) / error(h/2) ~ 4
    #
    # Use a somewhat loose threshold so the test isn't brittle.
    assert error_fine < error_coarse

    assert ratio > 3.0, (
        f"Expected approximately second-order L2 convergence; "
        f"got error ratio {ratio:.3f}"
    )


