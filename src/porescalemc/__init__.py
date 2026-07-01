"""
porescalemc — Multilevel Monte Carlo for pore-scale diffusion/permeability UQ.

The simplest way to run a study:

    from porescalemc import run_study, list_solvers

    list_solvers()                          # see what's available
    mean, err = run_study(
        solver="spectral_diffusion",
        n_grains=20, mu=0.08,               # geometry
        n_levels=3, m0=20, tolerance=0.05,  # MLMC
        resolution=16, n_directions=1,      # solver
    )

For more control, pass full config objects:

    from porescalemc import run_study
    from porescalemc.config import PackingConfig, MLMCConfig, SpectralConfig

    mean, err = run_study(
        solver="spectral_diffusion",
        packing_config=PackingConfig(n_grains=30, mu=0.07, periodic=True),
        mlmc_config=MLMCConfig(n_levels=4, m0=32, tolerance=0.02),
        solver_config=SpectralConfig(resolution=24, eta=1e-5, n_directions=3),
    )

Advanced: build factories yourself and use MLMCEstimator directly:

    from porescalemc import MLMCEstimator
    from porescalemc.config import PackingConfig, MLMCConfig, SpectralConfig, spectral_hierarchy_level
    from porescalemc.geometry.placement import sample_packing_from_name
    from porescalemc.solvers.spectral import SpectralDiffusionSolver

    def geo_factory(level, name, pcfg):
        return sample_packing_from_name(level, name, pcfg)

    def solver_factory(packing, level=0):
        cfg = spectral_hierarchy_level(SpectralConfig(resolution=12), level)
        return SpectralDiffusionSolver(spectral_config=cfg)

    est = MLMCEstimator(
        name="diffusion-uq",
        packing_config=PackingConfig(n_grains=20),
        mlmc_config=MLMCConfig(n_levels=3, m0=16),
        geometry_factory=geo_factory,
        solver_factory=solver_factory,
    )
    mean, err = est.run()
    print(est.diagnostics())
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from porescalemc.config import (
    FourierConfig,
    MLMCConfig,
    PackingConfig,
    SpectralConfig,
    spectral_hierarchy_level,
)
from porescalemc.geometry.fourier_field import (
    plot_field_slices,
    save_spectral_fields_to_vti,
    save_structured_vti,
    save_structured_vti_fields,
)
from porescalemc.geometry.placement import sample_packing_from_name
from porescalemc.mlmc.estimator import MLMCEstimator
from porescalemc.mlmc.statistics import (
    estimate_mlmc_work,
    estimate_observed_rates,
    mlmc_diagnostics_report,
    mlmc_speedup_report,
)
from porescalemc.solvers.fourier import FourierSolver
from porescalemc.solvers.packing import PackingStatsSolver
from porescalemc.solvers.registry import SOLVER_REGISTRY, available_qois, available_solvers, solver_class
from porescalemc.solvers.spectral import (
    SpectralAdvectionDiffusionSolver,
    SpectralDiffusionSolver,
    SpectralStokesSolver,
)
from porescalemc.solvers.tet import TetDiffusionSolver
from porescalemc.solvers.voxel import VoxelDiffusionSolver

try:
    from porescalemc.geometry.mesher import GmshPackingMesher
    from porescalemc.solvers.gmsh_mesher import GmshMesherQoISolver
except ImportError:
    pass

if TYPE_CHECKING:
    from numpy.typing import NDArray

__version__ = "0.1.0"
__author__ = "Matteo Icardi"


# ---------------------------------------------------------------------------
# User-friendly entry points
# ---------------------------------------------------------------------------

def list_solvers() -> dict[str, dict]:
    """Return a dict of solver metadata and print a human-readable table.

    Each value has keys: description, qoi_names, optional_dependency.

    Example
    -------
    >>> list_solvers()
    packing                  Algebraic geometry descriptors of the packing.
      QoIs: porosity, solid_fraction, ...
    spectral_diffusion       Pseudo-spectral VPM diffusion cell problem.
      QoIs: D_eff_x, D_eff_y, D_eff_z, ...
    ...
    """
    result: dict[str, dict] = {}
    for name, info in SOLVER_REGISTRY.items():
        result[name] = {
            "description": info.description,
            "qoi_names": info.qoi_names,
            "optional_dependency": info.optional_dependency,
        }
        req = f"  [requires: {info.optional_dependency}]" if info.optional_dependency else ""
        print(f"{name:<35} {info.description}{req}")
        print(f"  QoIs: {', '.join(info.qoi_names)}")
    return result


def list_qois(solver_name: str | None = None) -> dict[str, tuple]:
    """Return QoI names for one solver (or all solvers) and print them.

    Parameters
    ----------
    solver_name : str, optional
        Name of a specific solver. Omit to list QoIs for all solvers.

    Returns
    -------
    dict mapping solver name → tuple of QoI name strings.
    """
    if solver_name is not None:
        if solver_name not in SOLVER_REGISTRY:
            raise KeyError(
                f"Unknown solver '{solver_name}'. "
                f"Available: {list(SOLVER_REGISTRY)}"
            )
        names = SOLVER_REGISTRY[solver_name].qoi_names
        print(f"{solver_name}: {', '.join(names)}")
        return {solver_name: names}
    out: dict[str, tuple] = {}
    for name, info in SOLVER_REGISTRY.items():
        print(f"{name}: {', '.join(info.qoi_names)}")
        out[name] = info.qoi_names
    return out


def run_study(
    solver: str = "spectral_diffusion",
    *,
    # Geometry shortcuts (override with packing_config= for full control)
    n_grains: int = 20,
    mu: float = 0.08,
    coeff_var: float = 0.1,
    box_size: float = 1.0,
    psd: str = "lognormal",
    # MLMC shortcuts (override with mlmc_config= for full control)
    n_levels: int = 3,
    m0: int = 20,
    tolerance: float = 0.05,
    hierarchy: str = "g",
    n_workers: int = 1,
    algorithm: str = "fixed",
    # Solver shortcuts — meaningful for spectral/voxel solvers
    resolution: float = 16.0,
    n_directions: int = 1,
    eta: float = 1e-5,
    # Full config overrides
    packing_config: PackingConfig | None = None,
    mlmc_config: MLMCConfig | None = None,
    solver_config: object | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Run an MLMC study and return (mean, error_bound).

    Parameters
    ----------
    solver : str
        Solver name from ``list_solvers()``.  Defaults to
        ``"spectral_diffusion"``.
    n_grains : int
        Number of grains per unit-volume box.
    mu : float
        Mean grain radius.
    coeff_var : float
        Coefficient of variation of the grain radius PSD (0 = monodisperse).
    box_size : float
        Side length of the cubic periodic box.
    psd : str
        Grain size distribution: ``'lognormal'``, ``'uniform'``, or
        ``'constant'``.
    n_levels : int
        Number of MLMC levels (coarse to fine).
    m0 : int
        Samples at the coarsest level.
    tolerance : float
        Target relative RMSE for the MLMC estimator.
    hierarchy : str
        MLMC hierarchy type: ``'g'`` (grid refinement), ``'d'`` (domain
        scaling), ``'n'`` (grain count scaling), or combinations.
    n_workers : int
        Parallel workers (1 = serial).
    algorithm : str
        ``'fixed'`` or ``'adaptive'`` (Giles 2015).
    resolution : float
        FFT grid resolution (grid points per unit length per direction).
        Used by spectral and voxel solvers.
    n_directions : int
        Number of coordinate directions to solve (1–3).  3 returns the full
        3×3 effective tensor.
    eta : float
        VPM penalisation parameter for spectral solvers.
    packing_config : PackingConfig, optional
        Full packing configuration.  Overrides the geometry shortcuts.
    mlmc_config : MLMCConfig, optional
        Full MLMC configuration.  Overrides the MLMC shortcuts.
    solver_config : object, optional
        Solver-specific configuration (e.g. ``SpectralConfig``,
        ``FourierConfig``).  Overrides the solver shortcuts.

    Returns
    -------
    mean : ndarray  — estimated QoI mean vector
    error : ndarray — 99 % confidence half-width per QoI

    Examples
    --------
    Minimal (uses all defaults):

        mean, err = run_study()

    Spectral diffusion, 3 levels, full tensor:

        mean, err = run_study(
            solver="spectral_diffusion",
            n_grains=20, mu=0.08,
            n_levels=3, m0=20, tolerance=0.05,
            resolution=16, n_directions=3,
        )

    Voxel FV diffusion:

        mean, err = run_study(
            solver="voxel_diffusion",
            n_grains=25, mu=0.07,
            resolution=20, n_directions=1,
        )

    Stokes permeability:

        mean, err = run_study(
            solver="spectral_stokes",
            n_grains=15, mu=0.10,
            resolution=20, eta=1e-5,
        )
    """
    if solver not in SOLVER_REGISTRY:
        raise KeyError(
            f"Unknown solver '{solver}'. "
            f"Call list_solvers() to see available options."
        )

    # --- Build packing config ---
    if packing_config is None:
        packing_config = PackingConfig(
            n_grains=n_grains,
            mu=mu,
            coeff_var=coeff_var,
            xlen=box_size,
            ylen=box_size,
            zlen=box_size,
            psd=psd,
            periodic=True,
        )

    # --- Build MLMC config ---
    if mlmc_config is None:
        mlmc_config = MLMCConfig(
            n_levels=n_levels,
            m0=m0,
            tolerance=tolerance,
            hierarchy=hierarchy,
            n_workers=n_workers,
            algorithm=algorithm,
        )

    # --- Build solver factory ---
    solver_factory = _make_solver_factory(solver, solver_config, resolution, n_directions, eta)

    def geo_factory(level, name, pcfg):
        return sample_packing_from_name(level, name, pcfg)

    est = MLMCEstimator(
        name=f"porescalemc_{solver}",
        packing_config=packing_config,
        mlmc_config=mlmc_config,
        geometry_factory=geo_factory,
        solver_factory=solver_factory,
    )
    return est.run()


def _make_solver_factory(
    solver_name: str,
    solver_config: object | None,
    resolution: float,
    n_directions: int,
    eta: float,
):
    """Return a solver_factory function for the named solver."""
    from porescalemc.config import spectral_hierarchy_level

    if solver_name in ("spectral_diffusion", "spectral_stokes", "spectral_advection_diffusion"):
        base_cfg = solver_config or SpectralConfig(
            resolution=resolution,
            eta=eta,
            n_directions=n_directions,
        )
        cls = solver_class(solver_name)

        def factory(packing, level=0):
            cfg = spectral_hierarchy_level(base_cfg, level)
            return cls(spectral_config=cfg)

        return factory

    if solver_name == "voxel_diffusion":
        fcfg = solver_config or FourierConfig(resolution=resolution)

        def factory(packing, level=0):
            factor = 2.0 ** level
            lvl_cfg = FourierConfig(
                resolution=fcfg.resolution * factor,
                smoothing=fcfg.smoothing,
                smoothing_length=fcfg.smoothing_length,
            )
            return VoxelDiffusionSolver(
                fourier_config=lvl_cfg,
                n_directions=n_directions,
            )

        return factory

    if solver_name == "fourier":
        fcfg = solver_config or FourierConfig(resolution=resolution)

        def factory(packing, level=0):
            return FourierSolver(packing, fourier_config=fcfg)

        return factory

    if solver_name == "packing":
        def factory(packing, level=0):
            return PackingStatsSolver(packing)

        return factory

    if solver_name == "tet_diffusion":
        mesh_size = solver_config if isinstance(solver_config, float) else 0.05

        def factory(packing, level=0):
            size = mesh_size / (2.0 ** level)
            return TetDiffusionSolver(mesh_size=size, periodic=True)

        return factory

    # Generic fallback: instantiate solver class with no extra config
    cls = solver_class(solver_name)

    def factory(packing, level=0):
        return cls(packing)

    return factory


# ---------------------------------------------------------------------------
# Legacy convenience wrapper (kept for backwards compatibility)
# ---------------------------------------------------------------------------

def run_mlmc(
    packing_config: PackingConfig,
    mlmc_config: MLMCConfig,
    *,
    qoi_names: tuple[str, ...] | None = None,
    n_workers: int = 1,
    algorithm: str | None = None,
) -> tuple["np.ndarray", "np.ndarray"]:
    """High-level wrapper for the packing-statistics MLMC case.

    For a broader interface (choose solver, QoIs, etc.) use ``run_study()``.
    Returns (mean, error_bound).
    """
    if algorithm is not None:
        mlmc_config = MLMCConfig(**{**mlmc_config.__dict__, "algorithm": algorithm})
    if n_workers != 1:
        mlmc_config = MLMCConfig(**{**mlmc_config.__dict__, "n_workers": n_workers})

    qois = qoi_names or (
        "porosity", "solid_fraction", "specific_surface_area", "mean_equivalent_radius"
    )

    def geo(level, name, cfg):
        return sample_packing_from_name(level, name, cfg)

    def sol(packing):
        return PackingStatsSolver(packing, qoi_names=qois)

    est = MLMCEstimator(
        name="run_mlmc",
        packing_config=packing_config,
        mlmc_config=mlmc_config,
        geometry_factory=geo,
        solver_factory=sol,
        qoi_fn=lambda raw, p: raw,
    )
    return est.run()


def setup_logging(level: str = "INFO") -> None:
    """One-liner for basic logging during development or debugging."""
    import logging
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )
