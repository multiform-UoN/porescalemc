"""Voxel finite-volume diffusion solver.

This solver uses the Fourier reconstructed solid-fraction field as a voxel
geometry and solves the periodic scalar diffusion cell problem with a
matrix-free finite-volume operator.  It is slower and cruder than the spectral
VPM solver, but it provides an independent discretisation path for MLMC
experiments and convergence checks.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse.linalg import LinearOperator, cg

from porescalemc.config import FourierConfig
from porescalemc.geometry.fourier_field import porosity_field
from porescalemc.geometry.grains import Packing
from porescalemc.solvers.base import SolverProtocol


VOXEL_DIFFUSION_QOI_NAMES = (
    "diffusivity_x",
    "diffusivity_y",
    "diffusivity_z",
)


class VoxelDiffusionSolver(SolverProtocol):
    """Effective diffusivity from a periodic voxel finite-volume cell solve.

    Parameters
    ----------
    fourier_config : FourierConfig, optional
        Controls the voxel resolution and smoothing used to build the geometry.
    n_directions : int
        Number of coordinate directions to solve, from 1 to 3.
    solid_diffusivity : float
        Diffusivity assigned to solid voxels.  Use a small positive value for a
        nearly impermeable solid while keeping the operator well-conditioned.
    max_iter, tol : int, float
        Conjugate-gradient controls.
    """

    qoi_names = VOXEL_DIFFUSION_QOI_NAMES

    def __init__(
        self,
        packing: Packing | None = None,
        fourier_config: FourierConfig | None = None,
        n_directions: int = 1,
        solid_diffusivity: float = 1e-6,
        max_iter: int = 500,
        tol: float = 1e-7,
    ):
        if n_directions < 1 or n_directions > 3:
            raise ValueError("n_directions must be 1, 2 or 3.")
        self.packing = packing
        self.fourier_config = fourier_config or FourierConfig()
        self.n_directions = n_directions
        self.solid_diffusivity = float(solid_diffusivity)
        self.max_iter = int(max_iter)
        self.tol = float(tol)
        self._result: np.ndarray | None = None
        self._fields: dict[str, np.ndarray] = {}
        self._diagnostics: dict[str, float] = {}

    def setup(self, packing: Packing) -> None:
        self.packing = packing
        self._result = None
        self._fields = {}
        self._diagnostics = {}

    def solve(self) -> np.ndarray:
        if self.packing is None:
            raise RuntimeError("VoxelDiffusionSolver.setup must be called first.")

        phi_s = porosity_field(self.packing, self.fourier_config)
        phi_f = 1.0 - phi_s
        diffusivity = np.maximum(phi_f, self.solid_diffusivity)
        nx, ny, nz = phi_s.shape
        Lx, Ly, Lz = self.packing.box
        spacing = (Lx / nx, Ly / ny, Lz / nz)

        face = _face_diffusivities(diffusivity)
        operator = _build_operator(face, spacing, diffusivity.shape)

        values = np.zeros(self.n_directions, dtype=float)
        fields: dict[str, np.ndarray] = {
            "solid_fraction": phi_s.copy(),
            "porosity": phi_f.copy(),
            "voxel_diffusivity": diffusivity.copy(),
        }
        for axis in range(self.n_directions):
            rhs = _rhs_for_axis(face, spacing, axis)
            rhs = rhs - rhs.mean()
            chi = _solve_cg(operator, rhs.ravel(), self.max_iter, self.tol)
            chi = chi.reshape(diffusivity.shape)
            chi -= chi.mean()
            values[axis] = max(_effective_diffusivity(face, chi, spacing, axis), 0.0)
            fields[f"voxel_corrector_{'xyz'[axis]}"] = chi.copy()

        self._result = values
        self._fields = fields
        self._diagnostics = {
            "porosity_mean": float(phi_f.mean()),
            "solid_fraction_mean": float(phi_s.mean()),
            "diffusivity_mean": float(values.mean()),
            "diffusivity_min": float(values.min()),
            "diffusivity_max": float(values.max()),
        }
        for axis, value in enumerate(values):
            self._diagnostics[f"diffusivity_{'xyz'[axis]}"] = float(value)
        min_value = self._diagnostics["diffusivity_min"]
        self._diagnostics["diffusivity_anisotropy"] = (
            self._diagnostics["diffusivity_max"] / min_value if min_value > 0 else float("inf")
        )
        return values

    def solution_fields(self) -> dict[str, np.ndarray]:
        return dict(self._fields)

    def diagnostics(self) -> dict[str, float]:
        return dict(self._diagnostics)

    def close(self) -> None:
        return None


def _face_diffusivities(d: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Arithmetic face diffusivities in +x, +y and +z directions."""
    return (
        0.5 * (d + np.roll(d, -1, axis=0)),
        0.5 * (d + np.roll(d, -1, axis=1)),
        0.5 * (d + np.roll(d, -1, axis=2)),
    )


def _build_operator(
    face: tuple[np.ndarray, np.ndarray, np.ndarray],
    spacing: tuple[float, float, float],
    shape: tuple[int, int, int],
) -> LinearOperator:
    dx, dy, dz = spacing
    dpx, dpy, dpz = face
    dmx = np.roll(dpx, 1, axis=0)
    dmy = np.roll(dpy, 1, axis=1)
    dmz = np.roll(dpz, 1, axis=2)
    n = int(np.prod(shape))

    def matvec(flat: np.ndarray) -> np.ndarray:
        u = flat.reshape(shape)
        au = (
            dpx * (u - np.roll(u, -1, axis=0)) / dx**2
            + dmx * (u - np.roll(u, 1, axis=0)) / dx**2
            + dpy * (u - np.roll(u, -1, axis=1)) / dy**2
            + dmy * (u - np.roll(u, 1, axis=1)) / dy**2
            + dpz * (u - np.roll(u, -1, axis=2)) / dz**2
            + dmz * (u - np.roll(u, 1, axis=2)) / dz**2
        )
        au -= au.mean()
        return au.ravel()

    return LinearOperator((n, n), matvec=matvec, dtype=float)


def _rhs_for_axis(
    face: tuple[np.ndarray, np.ndarray, np.ndarray],
    spacing: tuple[float, float, float],
    axis: int,
) -> np.ndarray:
    dplus = face[axis]
    dminus = np.roll(dplus, 1, axis=axis)
    return (dplus - dminus) / spacing[axis]


def _solve_cg(operator: LinearOperator, rhs: np.ndarray, max_iter: int, tol: float) -> np.ndarray:
    if np.linalg.norm(rhs) < 1e-14:
        return np.zeros_like(rhs)
    try:
        sol, info = cg(operator, rhs, maxiter=max_iter, rtol=tol, atol=0.0)
    except TypeError:
        sol, info = cg(operator, rhs, maxiter=max_iter, tol=tol)
    if info < 0:
        raise RuntimeError(f"Voxel diffusion CG failed with info={info}.")
    return np.asarray(sol, dtype=float)


def _effective_diffusivity(
    face: tuple[np.ndarray, np.ndarray, np.ndarray],
    corrector: np.ndarray,
    spacing: tuple[float, float, float],
    axis: int,
) -> float:
    """Return mean periodic face flux in the imposed unit-gradient direction."""
    h = spacing[axis]
    dplus = face[axis]
    forward_grad = (np.roll(corrector, -1, axis=axis) - corrector) / h
    return float(np.mean(dplus * (1.0 + forward_grad)))
