"""OpenFOAM case preparation utilities for porescalemc.

Bridges between porescalemc packing/solver outputs and the OpenFOAM case
templates in ``templates/openfoam/``.  No OpenFOAM binary is invoked here;
this module only writes the input files that OpenFOAM needs.

Key outputs generated
---------------------
porosity.3d
    Space-separated table ``x y z K K K`` (isotropic permeability or
    diffusivity) consumed by the custom ``myDarcy`` / ``myLaplacian`` solvers.
    Format: first line ``nx ny nz``, then nx*ny*nz data rows.

constant/polyMesh/domainsize
    OpenFOAM #include dictionary fragment with domain geometry ($x1,$x2,…),
    grid resolution ($xgrid,$ygrid,$zgrid), scale factor, and BC type tokens.

system/tolerance
    Single line ``tol <value>;`` that overrides the solver tolerance.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Sequence

import numpy as np

from porescalemc.config import FourierConfig
from porescalemc.geometry.fourier_field import porosity_field
from porescalemc.geometry.grains import Packing


# ---------------------------------------------------------------------------
# Low-level writers
# ---------------------------------------------------------------------------

def write_porosity_3d(
    packing: Packing,
    path: str | Path,
    fourier_config: FourierConfig | None = None,
    f0: float = 1.0,
    f1: float = 0.0,
) -> None:
    """Write a ``porosity.3d`` file for the OpenFOAM porous-media solver.

    The file format is::

        nx ny nz
        x  y  z  val val val
        ...

    where *val* = f0 * phi_f + f1 * phi_s  (defaults: f0=1, f1=0 → porosity).
    For permeability pass f0 = K_fluid / K_ref, f1 = K_solid / K_ref.

    Coordinates are centred at the domain origin: x ∈ [-Lx/2, Lx/2].

    Parameters
    ----------
    packing : Packing
    path : str or Path
        Output file path (including filename, e.g. ``case/constant/porosity.3d``).
    fourier_config : FourierConfig, optional
        Controls FFT grid resolution / smoothing used to voxelise the packing.
    f0, f1 : float
        Linear combination weights: written_value = f0 * phi_f + f1 * phi_s.
    """
    cfg = fourier_config or FourierConfig(resolution=20.0)
    phi_s = porosity_field(packing, cfg)
    phi_f = 1.0 - phi_s
    nx, ny, nz = phi_s.shape
    Lx, Ly, Lz = packing.box

    field = f0 * phi_f + f1 * phi_s

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write(f"{nx} {ny} {nz}\n")
        xs = np.linspace(-Lx / 2, Lx / 2, nx)
        ys = np.linspace(-Ly / 2, Ly / 2, ny)
        zs = np.linspace(-Lz / 2, Lz / 2, nz)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    v = float(field[i, j, k])
                    fh.write(f"{xs[i]:.6g} {ys[j]:.6g} {zs[k]:.6g} {v:.6g} {v:.6g} {v:.6g}\n")


def write_domainsize(
    packing: Packing,
    path: str | Path,
    grid_res: float = 50.0,
    scale: float = 1.0,
    bcs: Sequence[str] = ("fixedValue", "fixedValue", "symmetry", "wall"),
    tol: float = 1e-4,
    n_procs: int = 1,
) -> None:
    """Write the ``domainsize`` dictionary fragment included by ``blockMeshDict``.

    Parameters
    ----------
    packing : Packing
    path : str or Path
        Output path, e.g. ``case/constant/polyMesh/domainsize``.
    grid_res : float
        Number of cells per unit length (cells = grid_res * box_length).
    scale : float
        ``scalegrid`` factor applied inside OpenFOAM's ``convertToMeters``.
    bcs : sequence of 4 str
        Boundary condition types for [inlet, outlet, lateral, pores].
        Supported: ``"fixedValue"``, ``"symmetry"``, ``"symmetryPlane"``,
        ``"empty"``, ``"cyclic"``, ``"wall"``.
    tol : float
        Solver tolerance written to ``system/tolerance``.
    n_procs : int
        Number of MPI processes (used to compute ndomainsx/y/z).
    """
    Lx, Ly, Lz = packing.box
    xgrid = max(1, int(grid_res * Lx))
    ygrid = max(1, int(grid_res * Ly))
    zgrid = max(1, int(grid_res * Lz))

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    bc_labels = ["inletbc", "outletbc", "lateralbc", "poresbc"]
    periodic_bcs = {"symmetry", "symmetryPlane", "empty", "cyclic", "cyclicAMI"}

    lines = [
        f"x1 {-0.5 * Lx:.6g};",
        f"y1 {-0.5 * Ly:.6g};",
        f"z1 {-0.5 * Lz:.6g};",
        f"x2 {+0.5 * Lx:.6g};",
        f"y2 {+0.5 * Ly:.6g};",
        f"z2 {+0.5 * Lz:.6g};",
        f"xgrid {xgrid};",
        f"ygrid {ygrid};",
        f"zgrid {zgrid};",
        f"scalegrid {scale:.6g};",
    ]
    for bc, lab in zip(bcs, bc_labels):
        if bc in periodic_bcs:
            lines.append(f"{lab}mesh {bc};")
        elif bc == "fixedValue":
            mesh_type = "patch" if "let" in lab else "wall"
            lines.append(f"{lab}mesh {mesh_type};")
        else:
            lines.append(f"{lab}mesh {bc};")
    lines.append("#inputMode merge")
    with path.open("w") as fh:
        fh.write("\n".join(lines) + "\n")


def write_tolerance(tol: float, path: str | Path) -> None:
    """Write ``system/tolerance`` with a single ``tol <value>;`` line."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write(f"tol {tol:.6g};\n#inputMode merge\n")


# ---------------------------------------------------------------------------
# High-level case preparation
# ---------------------------------------------------------------------------

def prepare_case(
    template: str | Path,
    output_dir: str | Path,
    packing: Packing,
    fourier_config: FourierConfig | None = None,
    grid_res: float = 50.0,
    scale: float = 1.0,
    bcs: Sequence[str] = ("fixedValue", "fixedValue", "symmetry", "wall"),
    tol: float = 1e-4,
    n_procs: int = 1,
    overwrite: bool = False,
) -> Path:
    """Populate an OpenFOAM case directory from a template.

    Copies the template directory tree to *output_dir*, then writes the
    packing-derived files (``domainsize``, ``porosity.3d``, ``tolerance``).

    Parameters
    ----------
    template : str or Path
        Path to the template directory (e.g. ``templates/openfoam/darcy``).
    output_dir : str or Path
        Destination case directory.  Created if it does not exist.
    packing : Packing
    fourier_config : FourierConfig, optional
    grid_res : float
    scale : float
    bcs : sequence of 4 str
    tol : float
    n_procs : int
    overwrite : bool
        If True, remove *output_dir* before copying. Default False.

    Returns
    -------
    Path
        Absolute path to the prepared case directory.
    """
    template = Path(template).resolve()
    output_dir = Path(output_dir).resolve()

    if not template.is_dir():
        raise FileNotFoundError(f"Template not found: {template}")

    if overwrite and output_dir.exists():
        shutil.rmtree(output_dir)

    shutil.copytree(template, output_dir, dirs_exist_ok=True)

    write_domainsize(
        packing,
        output_dir / "constant" / "polyMesh" / "domainsize",
        grid_res=grid_res,
        scale=scale,
        bcs=bcs,
        tol=tol,
        n_procs=n_procs,
    )
    write_porosity_3d(
        packing,
        output_dir / "constant" / "porosity.3d",
        fourier_config=fourier_config,
    )
    write_tolerance(tol, output_dir / "system" / "tolerance")

    return output_dir


# ---------------------------------------------------------------------------
# Template introspection
# ---------------------------------------------------------------------------

_REQUIRED_FILES = [
    "system/controlDict",
    "system/fvSchemes",
    "system/fvSolution",
    "constant",
]


def validate_template(template: str | Path) -> list[str]:
    """Return a list of missing required paths in a template directory.

    Empty list means the template is structurally complete.
    """
    template = Path(template)
    missing = []
    for rel in _REQUIRED_FILES:
        if not (template / rel).exists():
            missing.append(rel)
    return missing
