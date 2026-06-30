"""User-facing registry of built-in solver adapters and QoI names."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from porescalemc.solvers.fourier import FourierSolver
from porescalemc.solvers.packing import PACKING_QOI_NAMES, PackingStatsSolver
from porescalemc.solvers.spectral import (
    SpectralAdvectionDiffusionSolver,
    SpectralDiffusionSolver,
    SpectralStokesSolver,
)
from porescalemc.solvers.tet import TET_DIFFUSION_QOI_NAMES, TetDiffusionSolver
from porescalemc.solvers.voxel import VOXEL_DIFFUSION_QOI_NAMES, VoxelDiffusionSolver


@dataclass(frozen=True)
class SolverInfo:
    """Metadata for a built-in solver adapter."""

    name: str
    cls: type
    qoi_names: tuple[str, ...]
    description: str
    optional_dependency: str | None = None


SOLVER_REGISTRY: dict[str, SolverInfo] = {
    "packing": SolverInfo(
        name="packing",
        cls=PackingStatsSolver,
        qoi_names=PACKING_QOI_NAMES,
        description="Algebraic geometry descriptors of the packing.",
    ),
    "fourier": SolverInfo(
        name="fourier",
        cls=FourierSolver,
        qoi_names=("mean_solid_fraction", "upscaled_permeability"),
        description="Fourier-field porosity and permeability closure QoIs.",
    ),
    "spectral_diffusion": SolverInfo(
        name="spectral_diffusion",
        cls=SpectralDiffusionSolver,
        qoi_names=("D_eff_x", "D_eff_y", "D_eff_z", "D_tensor_3x3_if_n_directions_3"),
        description="Pseudo-spectral VPM diffusion cell problem.",
    ),
    "spectral_stokes": SolverInfo(
        name="spectral_stokes",
        cls=SpectralStokesSolver,
        qoi_names=("K_eff_x", "K_eff_y", "K_eff_z", "K_tensor_3x3_if_n_directions_3"),
        description="Pseudo-spectral Brinkman-Stokes permeability problem.",
    ),
    "spectral_advection_diffusion": SolverInfo(
        name="spectral_advection_diffusion",
        cls=SpectralAdvectionDiffusionSolver,
        qoi_names=("D_adv_x", "D_adv_y", "D_adv_z", "K_eff_x", "K_eff_y", "K_eff_z"),
        description="Pseudo-spectral finite-Pe advection-diffusion estimate.",
    ),
    "voxel_diffusion": SolverInfo(
        name="voxel_diffusion",
        cls=VoxelDiffusionSolver,
        qoi_names=VOXEL_DIFFUSION_QOI_NAMES,
        description="Periodic voxel finite-volume diffusion cell problem.",
    ),
    "tet_diffusion": SolverInfo(
        name="tet_diffusion",
        cls=TetDiffusionSolver,
        qoi_names=TET_DIFFUSION_QOI_NAMES,
        description="Gmsh tetrahedral mesh-backed Bruggeman diffusivity estimate.",
        optional_dependency="gmsh",
    ),
}


def available_solvers() -> tuple[str, ...]:
    """Return names of built-in solver adapters."""
    return tuple(SOLVER_REGISTRY)


def available_qois(solver_name: str | None = None) -> dict[str, tuple[str, ...]] | tuple[str, ...]:
    """Return QoI names for one solver or all solvers."""
    if solver_name is not None:
        return SOLVER_REGISTRY[solver_name].qoi_names
    return {name: info.qoi_names for name, info in SOLVER_REGISTRY.items()}


def solver_class(name: str) -> type:
    """Return the class implementing a built-in solver adapter."""
    return SOLVER_REGISTRY[name].cls
