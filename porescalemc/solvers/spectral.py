"""Spectral VPM solvers for periodic porous media.

Three solvers implement SolverProtocol:

SpectralDiffusionSolver
    Effective diffusivity tensor D_eff/D0 via the periodic VPM cell problem.
    Uses preconditioned CG (operator is SPD). Supports Dirichlet (absorbing
    solid) and Neumann (no-flux solid) immersed boundary conditions.

SpectralStokesSolver
    Effective permeability tensor K_eff via VPM Brinkman-Stokes Richardson
    iteration with solenoidal projection. Improved damping parameter lambda
    prevents blow-up near the Nyquist frequency.

SpectralAdvectionDiffusionSolver
    Effective dispersion tensor for transport in Stokes flow at finite Peclet
    number. Uses outer Picard iteration with CG inner solves.

Mathematical background
-----------------------
All solvers use the Volume Penalisation Method (VPM). A solid indicator
field phi_s is reconstructed from the grain packing via analytical Fourier
transforms (see geometry/fourier_field.py). The fluid field is phi_f = 1 - phi_s.

Diffusion cell problem (Dirichlet VPM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Decompose c(x) = xi/Li + c_tilde(x) where xi is the coordinate along
direction i, c_tilde is a periodic fluctuation, and 1/Li is the imposed
mean concentration gradient. The VPM equation -D0∇²c + (phi_s/η)c = 0 gives:

    A(c_tilde) = -(phi_s/η) * xi/Li,     A = -D0∇² + (phi_s/η)

A is SPD → CG with Fourier preconditioner. Effective diffusivity:

    D_eff_i/D0 = phi_f_mean + Li * mean(phi_f * ∂c_tilde/∂xi)

Empty domain limit: phi_s = 0 → RHS = 0 → c_tilde = 0 → D_eff = D0. ✓

Diffusion cell problem (Neumann VPM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Modified diffusivity D(x) = D0 * phi_f(x). The cell problem becomes:

    A_N(c_tilde) = -D0/Li * ∂phi_s/∂xi,   A_N = -∇·(phi_f D0 ∇·) + ε_reg

A_N is symmetric positive semidefinite; ε_reg regularises the zero eigenvalue.

Stokes (Richardson + solenoidal projection)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Drive with body force f = e_i (unit vector in direction i), measure:
    K_eff_i = mean(u_i)   (assumes μ = 1)

Richardson iteration:
    u^{n+1} = P_sol[ (λ u^n - (phi_s/η) u^n + f) / (μ 4π²k² + λ) ]

where P_sol is the solenoidal (divergence-free) projector in k-space and
lambda is chosen to exceed the maximum eigenvalue of phi_s/η.
"""

from __future__ import annotations

import logging
from typing import Callable

import numpy as np

from porescalemc.config import FourierConfig, SpectralConfig
from porescalemc.geometry.fourier_field import porosity_field
from porescalemc.geometry.grains import Packing
from porescalemc.solvers.base import SolverProtocol

_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# FFT helpers
# ---------------------------------------------------------------------------

def _build_k_arrays(packing: Packing, nx: int, ny: int, nz: int):
    """Wavenumber, grid-spacing, and coordinate arrays on the FFT grid.

    Returns
    -------
    kx, ky, kz : ndarray
        Broadcastable wavenumber arrays, shapes (nx,1,1), (1,ny,1), (1,1,nz).
    k2 : ndarray, shape (nx, ny, nz)
        |k|^2.
    dx, dy, dz : float
        Grid spacings.
    xi, yi, zi : ndarray
        Physical coordinate arrays (0 to L-dL), shapes (nx,1,1) etc.
        Not periodic; used only for constructing the cell-problem RHS.
    """
    Lx, Ly, Lz = packing.box
    dx = Lx / nx
    dy = Ly / ny
    dz = Lz / nz
    kx = np.fft.fftfreq(nx, d=dx).reshape(nx, 1, 1)
    ky = np.fft.fftfreq(ny, d=dy).reshape(1, ny, 1)
    kz = np.fft.fftfreq(nz, d=dz).reshape(1, 1, nz)
    k2 = kx ** 2 + ky ** 2 + kz ** 2
    xi = np.arange(nx).reshape(nx, 1, 1) * dx
    yi = np.arange(ny).reshape(1, ny, 1) * dy
    zi = np.arange(nz).reshape(1, 1, nz) * dz
    return kx, ky, kz, k2, dx, dy, dz, xi, yi, zi


def _laplacian_k(c: np.ndarray, k2: np.ndarray) -> np.ndarray:
    """Spectral Laplacian: ∇²c = IFFT(-4π²|k|² FFT(c))."""
    return np.real(np.fft.ifftn(-4.0 * np.pi ** 2 * k2 * np.fft.fftn(c)))


def _gradient_k(c: np.ndarray, kx, ky, kz):
    """Spectral gradient: (∂c/∂x, ∂c/∂y, ∂c/∂z) = IFFT(2πi k FFT(c))."""
    c_hat = np.fft.fftn(c)
    tpi = 2.0j * np.pi
    gx = np.real(np.fft.ifftn(tpi * kx * c_hat))
    gy = np.real(np.fft.ifftn(tpi * ky * c_hat))
    gz = np.real(np.fft.ifftn(tpi * kz * c_hat))
    return gx, gy, gz


def _divergence_k(vx, vy, vz, kx, ky, kz) -> np.ndarray:
    """Spectral divergence: ∇·v = IFFT(2πi (kx V̂x + ky V̂y + kz V̂z))."""
    tpi = 2.0j * np.pi
    return np.real(np.fft.ifftn(
        tpi * kx * np.fft.fftn(vx) +
        tpi * ky * np.fft.fftn(vy) +
        tpi * kz * np.fft.fftn(vz)
    ))


# ---------------------------------------------------------------------------
# Preconditioned CG
# ---------------------------------------------------------------------------

def _cg_solve(
    A_op: Callable,
    b: np.ndarray,
    x0: np.ndarray | None = None,
    max_iter: int = 500,
    tol: float = 1e-6,
    M_op: Callable | None = None,
) -> np.ndarray:
    """Preconditioned Conjugate Gradient for A x = b, A symmetric positive definite.

    A_op and M_op operate on flat 1-D float arrays and return flat 1-D arrays.
    M_op is the preconditioner inverse application: z = M⁻¹ r.
    Converges in at most min(n, max_iter) steps for exact arithmetic.
    """
    x = np.zeros(b.size) if x0 is None else x0.copy()
    r = b - A_op(x)
    z = M_op(r) if M_op is not None else r.copy()
    p = z.copy()
    rz = float(np.dot(r, z))
    b_norm = float(np.linalg.norm(b)) + 1e-15

    for _ in range(max_iter):
        Ap = A_op(p)
        pAp = float(np.dot(p, Ap))
        if abs(pAp) < 1e-30:
            break
        alpha = rz / pAp
        x += alpha * p
        r -= alpha * Ap
        if float(np.linalg.norm(r)) / b_norm < tol:
            break
        z_new = M_op(r) if M_op is not None else r.copy()
        rz_new = float(np.dot(r, z_new))
        beta = rz_new / max(abs(rz), 1e-30)
        p = z_new + beta * p
        rz = rz_new

    return x


# ---------------------------------------------------------------------------
# Diffusion solver
# ---------------------------------------------------------------------------

class SpectralDiffusionSolver(SolverProtocol):
    """Effective diffusivity tensor D_eff/D0 via periodic VPM cell problem.

    For each coordinate direction i, solves:

    Dirichlet VPM (bc_solid='dirichlet'):
        A(c_tilde) = -(phi_s/η) * xi/Li,    A = -D0∇² + (phi_s/η)

    Neumann VPM (bc_solid='neumann'):
        A_N(c_tilde) = -D0/Li * ∂phi_s/∂xi,  A_N = -∇·(phi_f D0 ∇·) + ε_reg

    Both operators are SPD and solved by CG with a Fourier preconditioner.

    Parameters
    ----------
    spectral_config : SpectralConfig, optional
        Solver parameters. If None, uses SpectralConfig() defaults.
    fourier_config : FourierConfig, optional
        Porosity field parameters. If None, derived from spectral_config.
    packing : Packing, optional
        Can be supplied here or via setup().
    eta, max_iter, tol : legacy kwargs
        Override the corresponding spectral_config fields (backward compat).

    Returns (solve)
    ---------------
    np.ndarray, shape (n_directions,)
        D_eff_i / D0 for i = 0 … n_directions - 1.
    """

    def __init__(
        self,
        spectral_config: SpectralConfig | None = None,
        fourier_config: FourierConfig | None = None,
        packing: Packing | None = None,
    ):
        self.cfg = spectral_config or SpectralConfig()
        self.fourier_config = fourier_config or FourierConfig(
            resolution=self.cfg.resolution,
            smoothing=self.cfg.smoothing,
            smoothing_length=self.cfg.smoothing_length,
        )
        self.packing: Packing | None = packing
        self._result: np.ndarray | None = None

    def setup(self, packing: Packing) -> None:
        self.packing = packing
        self._result = None

    def solve(self) -> np.ndarray:
        if self.packing is None:
            raise RuntimeError("setup(packing) must be called first.")

        phi_s = porosity_field(self.packing, self.fourier_config)
        phi_f = 1.0 - phi_s
        nx, ny, nz = phi_s.shape

        kx, ky, kz, k2, dx, dy, dz, xi, yi, zi = _build_k_arrays(
            self.packing, nx, ny, nz
        )
        Lx, Ly, Lz = self.packing.box
        phi_f_mean = float(phi_f.mean())

        coords = [(xi, kx, Lx), (yi, ky, Ly), (zi, kz, Lz)]
        D_eff = np.zeros(self.cfg.n_directions)

        for i, (coord, ki, Li) in enumerate(coords[: self.cfg.n_directions]):
            if self.cfg.bc_solid == "neumann":
                D_eff[i] = self._solve_neumann(
                    phi_s, phi_f, phi_f_mean, coord, ki, Li,
                    kx, ky, kz, k2, nx, ny, nz,
                )
            else:
                D_eff[i] = self._solve_dirichlet(
                    phi_s, phi_f, phi_f_mean, coord, ki, Li, k2, nx, ny, nz,
                )

        self._result = D_eff
        return D_eff

    # ------------------------------------------------------------------
    # Internal solves

    def _solve_dirichlet(
        self, phi_s, phi_f, phi_f_mean,
        coord, ki, Li, k2, nx, ny, nz,
    ) -> float:
        eta = self.cfg.eta
        pen = phi_s / eta
        # Fourier preconditioner: approximate A with constant-coefficient operator
        lam_pc = max(float(phi_s.mean()) / eta, 1.0 / eta * 1e-3, 1e-6)

        def A_op(c_flat: np.ndarray) -> np.ndarray:
            c3d = c_flat.reshape(nx, ny, nz)
            return (-_laplacian_k(c3d, k2) + pen * c3d).ravel()

        def M_op(r_flat: np.ndarray) -> np.ndarray:
            r3d = r_flat.reshape(nx, ny, nz)
            return np.real(
                np.fft.ifftn(np.fft.fftn(r3d) / (4.0 * np.pi ** 2 * k2 + lam_pc))
            ).ravel()

        b = (-pen * coord / Li).ravel()
        c_tilde = _cg_solve(
            A_op, b, max_iter=self.cfg.max_iter, tol=self.cfg.tol, M_op=M_op,
        )
        c3d = c_tilde.reshape(nx, ny, nz)
        dc = np.real(np.fft.ifftn(2.0j * np.pi * ki * np.fft.fftn(c3d)))
        return max(phi_f_mean + Li * float(np.mean(phi_f * dc)), 0.0)

    def _solve_neumann(
        self, phi_s, phi_f, phi_f_mean,
        coord, ki, Li, kx, ky, kz, k2, nx, ny, nz,
    ) -> float:
        eps_reg = 1e-8
        phi_f_avg = max(phi_f_mean, 1e-6)

        def A_op(c_flat: np.ndarray) -> np.ndarray:
            c3d = c_flat.reshape(nx, ny, nz)
            gx, gy, gz = _gradient_k(c3d, kx, ky, kz)
            Ac = -_divergence_k(
                phi_f * gx, phi_f * gy, phi_f * gz, kx, ky, kz
            )
            return (Ac + eps_reg * c3d).ravel()

        def M_op(r_flat: np.ndarray) -> np.ndarray:
            r3d = r_flat.reshape(nx, ny, nz)
            denom = phi_f_avg * 4.0 * np.pi ** 2 * k2 + eps_reg
            result = np.real(np.fft.ifftn(np.fft.fftn(r3d) / denom))
            result -= result.mean()
            return result.ravel()

        # RHS: -1/Li * ∂phi_s/∂xi (spectral gradient of phi_s)
        dphi_s = np.real(np.fft.ifftn(2.0j * np.pi * ki * np.fft.fftn(phi_s)))
        b = (-dphi_s / Li).ravel()
        c_tilde = _cg_solve(
            A_op, b, max_iter=self.cfg.max_iter, tol=self.cfg.tol, M_op=M_op,
        )
        c3d = c_tilde.reshape(nx, ny, nz)
        dc = np.real(np.fft.ifftn(2.0j * np.pi * ki * np.fft.fftn(c3d)))
        return max(phi_f_mean + Li * float(np.mean(phi_f * dc)), 0.0)

    def close(self) -> None:
        pass


# ---------------------------------------------------------------------------
# Stokes solver
# ---------------------------------------------------------------------------

class SpectralStokesSolver(SolverProtocol):
    """Effective permeability tensor K_eff via VPM Brinkman-Stokes (Richardson).

    For each direction i, drives flow with constant body force f = e_i and
    measures K_eff_i = mean(u_i), assuming viscosity μ = 1.

    The damping parameter λ is chosen to exceed the maximum eigenvalue of the
    penalisation operator phi_s/η and the viscous term at the Nyquist frequency:

        λ = max(phi_s_max / η,  μ (2π k_Nyq)² / 4) with floor 1.0

    Parameters
    ----------
    spectral_config : SpectralConfig, optional
    fourier_config : FourierConfig, optional
    packing : Packing, optional
    mu : float
        Dynamic viscosity (default 1.0).
    eta, max_iter, tol : legacy kwargs

    Returns (solve)
    ---------------
    np.ndarray, shape (n_directions,)
        K_eff_i for i = 0 … n_directions - 1.
    """

    def __init__(
        self,
        spectral_config: SpectralConfig | None = None,
        fourier_config: FourierConfig | None = None,
        packing: Packing | None = None,
        mu: float = 1.0,
    ):
        self.cfg = spectral_config or SpectralConfig()
        self.mu = mu
        self.fourier_config = fourier_config or FourierConfig(
            resolution=self.cfg.resolution,
            smoothing=self.cfg.smoothing,
            smoothing_length=self.cfg.smoothing_length,
        )
        self.packing: Packing | None = packing
        self._result: np.ndarray | None = None
        self._velocity_fields: list[np.ndarray] | None = None  # u per direction

    def setup(self, packing: Packing) -> None:
        self.packing = packing
        self._result = None
        self._velocity_fields = None

    def solve(self) -> np.ndarray:
        if self.packing is None:
            raise RuntimeError("setup(packing) must be called first.")

        phi_s = porosity_field(self.packing, self.fourier_config)
        nx, ny, nz = phi_s.shape
        N = nx * ny * nz
        Lx, Ly, Lz = self.packing.box

        kx, ky, kz, k2, dx, dy, dz, _, _, _ = _build_k_arrays(
            self.packing, nx, ny, nz
        )
        k2_safe = np.where(k2 > 1e-15, k2, 1.0)
        k_vec = np.stack(np.broadcast_arrays(kx, ky, kz), axis=0)  # (3,nx,ny,nz)

        def project_solenoidal(F_hat: np.ndarray) -> np.ndarray:
            """Remove compressible part: P = I - k⊗k/|k|²."""
            k_dot_F = (k_vec * F_hat).sum(axis=0, keepdims=True)
            return F_hat - k_vec * (k_dot_F / k2_safe)

        mu = self.mu
        eta = self.cfg.eta
        k_nyq = 0.5 / min(dx, dy, dz)
        lam = max(
            float(phi_s.max()) / eta,
            mu * (2.0 * np.pi * k_nyq) ** 2 * 0.25,
            1.0,
        )

        phi_pen = phi_s / eta
        denom = mu * 4.0 * np.pi ** 2 * k2 + lam  # shape (nx,ny,nz)

        K_eff = np.zeros(self.cfg.n_directions)
        velocity_fields: list[np.ndarray] = []

        for i in range(self.cfg.n_directions):
            rhs_hat = np.zeros((3, nx, ny, nz), dtype=complex)
            rhs_hat[i, 0, 0, 0] = float(N)  # unit body force in direction i

            u = np.zeros((3, nx, ny, nz), dtype=float)
            for _it in range(self.cfg.max_iter):
                term = (lam - phi_pen) * u
                term_hat = np.stack(
                    [np.fft.fftn(term[d]) for d in range(3)], axis=0
                )
                u_hat_new = project_solenoidal((term_hat + rhs_hat) / denom)
                u_new = np.stack(
                    [np.real(np.fft.ifftn(u_hat_new[d])) for d in range(3)], axis=0
                )
                diff = float(
                    np.linalg.norm(u_new - u) / (np.linalg.norm(u_new) + 1e-15)
                )
                u = u_new
                if diff < self.cfg.tol:
                    break

            K_eff[i] = float(u[i].mean())
            velocity_fields.append(u.copy())

        self._velocity_fields = velocity_fields
        self._result = K_eff
        return K_eff

    def close(self) -> None:
        self._velocity_fields = None


# ---------------------------------------------------------------------------
# Advection-diffusion solver
# ---------------------------------------------------------------------------

class SpectralAdvectionDiffusionSolver(SolverProtocol):
    """Effective dispersion tensor via VPM advection-diffusion cell problem.

    Combines a Stokes velocity field (computed internally) with the diffusion
    cell problem at finite Peclet number Pe = cfg.peclet. For Pe = 0, returns
    the same D_eff as SpectralDiffusionSolver (Dirichlet VPM).

    Algorithm
    ---------
    For each direction i:
    1. Solve Stokes with body force e_i → velocity field u(x).
    2. Solve outer Picard iteration treating Pe*u·∇c_tilde as a source:
           A(c_tilde^{k+1}) = b_diff - Pe*u·∇c_tilde^k - Pe*u_i/Li
       where A = -D0∇² + (phi_s/η) and b_diff = -(phi_s/η)*xi/Li.
    3. D_eff_i = phi_f_mean + Li * mean(phi_f * ∂c_tilde/∂xi).

    The advection term Pe*u·∇c renders the full system non-symmetric; the
    outer fixed-point splits it into a symmetric CG step (inner) + explicit
    advection correction (outer), which converges for Pe not too large.

    Returns (solve)
    ---------------
    np.ndarray, shape (2 * n_directions,)
        Concatenation of [D_eff_0, ..., D_eff_{n-1}, K_eff_0, ..., K_eff_{n-1}].
        K_eff values are the Stokes permeabilities used to generate the velocity.
    """

    def __init__(
        self,
        spectral_config: SpectralConfig | None = None,
        fourier_config: FourierConfig | None = None,
        packing: Packing | None = None,
    ):
        self.cfg = spectral_config or SpectralConfig()
        self.fourier_config = fourier_config or FourierConfig(
            resolution=self.cfg.resolution,
            smoothing=self.cfg.smoothing,
            smoothing_length=self.cfg.smoothing_length,
        )
        self.packing: Packing | None = packing
        self._result: np.ndarray | None = None

    def setup(self, packing: Packing) -> None:
        self.packing = packing
        self._result = None

    def solve(self) -> np.ndarray:
        if self.packing is None:
            raise RuntimeError("setup(packing) must be called first.")

        # Step 1: Stokes velocity fields
        stokes = SpectralStokesSolver(
            spectral_config=self.cfg,
            fourier_config=self.fourier_config,
        )
        stokes.setup(self.packing)
        K_eff = stokes.solve()
        vel_fields = stokes._velocity_fields  # list of (3,nx,ny,nz) arrays

        phi_s = porosity_field(self.packing, self.fourier_config)
        phi_f = 1.0 - phi_s
        nx, ny, nz = phi_s.shape
        N = nx * ny * nz
        Lx, Ly, Lz = self.packing.box

        kx, ky, kz, k2, dx, dy, dz, xi, yi, zi = _build_k_arrays(
            self.packing, nx, ny, nz
        )

        eta = self.cfg.eta
        Pe = self.cfg.peclet
        phi_f_mean = float(phi_f.mean())
        pen = phi_s / eta
        lam_pc = max(float(phi_s.mean()) / eta, 1e-6)

        def A_diff(c_flat: np.ndarray) -> np.ndarray:
            c3d = c_flat.reshape(nx, ny, nz)
            return (-_laplacian_k(c3d, k2) + pen * c3d).ravel()

        def M_pc(r_flat: np.ndarray) -> np.ndarray:
            r3d = r_flat.reshape(nx, ny, nz)
            return np.real(
                np.fft.ifftn(np.fft.fftn(r3d) / (4.0 * np.pi ** 2 * k2 + lam_pc))
            ).ravel()

        coords = [(xi, kx, Lx), (yi, ky, Ly), (zi, kz, Lz)]
        D_eff_ad = np.zeros(self.cfg.n_directions)

        for i, (coord, ki, Li) in enumerate(coords[: self.cfg.n_directions]):
            u_field = vel_fields[i]  # velocity field driven by e_i
            b_diff = (-pen * coord / Li).ravel()

            c_tilde = np.zeros(N)
            for _outer in range(50):
                c3d = c_tilde.reshape(nx, ny, nz)
                # Advection of fluctuation and of mean gradient field
                dc_dx, dc_dy, dc_dz = _gradient_k(c3d, kx, ky, kz)
                adv = Pe * (
                    u_field[0] * dc_dx +
                    u_field[1] * dc_dy +
                    u_field[2] * dc_dz
                )
                adv_mean = Pe * u_field[i] / Li
                b_total = b_diff - adv.ravel() - adv_mean.ravel()
                c_new = _cg_solve(
                    A_diff, b_total, x0=c_tilde,
                    max_iter=self.cfg.max_iter, tol=self.cfg.tol,
                    M_op=M_pc,
                )
                diff = float(
                    np.linalg.norm(c_new - c_tilde) / (np.linalg.norm(c_new) + 1e-15)
                )
                c_tilde = c_new
                if diff < self.cfg.tol * 10:
                    break

            c3d_f = c_tilde.reshape(nx, ny, nz)
            dc = np.real(np.fft.ifftn(2.0j * np.pi * ki * np.fft.fftn(c3d_f)))
            D_eff_ad[i] = max(phi_f_mean + Li * float(np.mean(phi_f * dc)), 0.0)

        self._result = np.concatenate([D_eff_ad, K_eff])
        return self._result

    def close(self) -> None:
        pass
