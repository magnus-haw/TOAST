# examples/manufactured.py

import numpy as np
import ufl

from mpi4py import MPI
from petsc4py import PETSc

from dolfinx import fem, mesh
from dolfinx.fem.petsc import LinearProblem

from toast.verification import (
    manufactured_temperature_ufl,
    manufactured_source_ufl,
)

from toast.visualization import plot_temperature_field


LX = 1.0
LY = 0.7

K = 0.15

T0 = 293.15
AMPLITUDE = 10.0


def solve_manufactured_problem(nx=40, ny=28):

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
        K
        * ufl.dot(
            ufl.grad(T),
            ufl.grad(v),
        )
        * ufl.dx
    )

    L = Q * v * ufl.dx

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
    T_h.name = "Temperature"

    problem = LinearProblem(
        a,
        L,
        bcs=[bc],
        u=T_h,
        petsc_options_prefix="manufactured_example_",
        petsc_options={
            "ksp_type": "preonly",
            "pc_type": "lu",
        },
    )

    problem.solve()
    T_h.x.scatter_forward()

    exact = manufactured_temperature_ufl(
        domain,
        lx=LX,
        ly=LY,
        t0=T0,
        amplitude=AMPLITUDE,
    )

    error_form = fem.form(
        (T_h - exact) ** 2 * ufl.dx
    )

    local_error_sq = fem.assemble_scalar(
        error_form
    )

    global_error_sq = domain.comm.allreduce(
        local_error_sq,
        op=MPI.SUM,
    )

    return domain, T_h, np.sqrt(global_error_sq)


if __name__ == "__main__":

    domain, temperature, error = (
        solve_manufactured_problem()
    )

    if domain.comm.rank == 0:
        print(f"L2 error = {error:.6e}")

    plot_temperature_field(
        domain,
        temperature,
        title="Manufactured 2-D solution",
    )
