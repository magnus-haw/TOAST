# Installation

## Recommended environment

The project is developed with Python 3.12 and the conda-forge DOLFINx stack. The repository includes `environment.yml`:

```yaml
name: toast

channels:
  - conda-forge

dependencies:
  - python=3.12
  - fenics-dolfinx
  - mpich
  - gmsh
  - python-gmsh
  - scipy
  - pandas
  - matplotlib
  - pyvista
  - ipykernel
  - pip
  - pytest
```

Create and activate the environment:

```bash
conda env create -f environment.yml
conda activate toast
```

Install the local package in editable mode:

```bash
pip install -e .
```

Editable installation means changes under `src/toast/` are available immediately without reinstalling after each edit.

## Verify the numerical stack

A useful low-level check is:

```bash
python - <<'PY'
import dolfinx
import gmsh
import mpi4py
import petsc4py
from petsc4py import PETSc

print("DOLFINx:", dolfinx.__version__)
print("Gmsh:   ", gmsh.__version__)
print("mpi4py: ", mpi4py.__version__)
print("petsc4py:", petsc4py.__version__)
print("PETSc:  ", PETSc.Sys.getVersion())
PY
```

Then run the project tests:

```bash
pytest -q
```

## Why Conda owns the compiled stack

`pyproject.toml` intentionally does not ask pip to install DOLFINx, PETSc, MPI, or Gmsh. These packages contain compiled libraries that must be mutually compatible. Installing them together from conda-forge reduces the chance of mixing incompatible MPI, PETSc, OpenMP, HDF5, or sparse-solver libraries.

## macOS and Homebrew OpenMP

If PETSc fails during import with an error involving a Homebrew `libomp.dylib`, inspect:

```bash
echo "$DYLD_LIBRARY_PATH"
```

A globally exported Homebrew OpenMP path can override the Conda environment's compatible runtime. For the `toast` environment, use:

```bash
unset DYLD_LIBRARY_PATH
```

and remove any global `export DYLD_LIBRARY_PATH=...homebrew...libomp...` from shell startup files unless another application genuinely requires it.

## Developer tools

The `pyproject.toml` defines optional Python development dependencies. If desired:

```bash
pip install -e ".[dev,visualization,analysis]"
```

This adds tools such as pytest, Ruff, mypy, build, PyVista, and pandas through pip. The compiled FEM stack should still remain Conda-managed.

Typical commands are:

```bash
pytest -q
ruff check src tests examples
ruff format --check src tests examples
```
