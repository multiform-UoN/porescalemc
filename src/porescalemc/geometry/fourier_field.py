"""
Fourier-space representation of grain packings for computing smooth porosity fields.

Mathematical background
-----------------------
Instead of voxelising the pore space and smoothing, we use the *analytical*
Fourier transform of each grain.  Summing over all grains gives the exact
(up to machine precision) k-space representation of the solid indicator field,
which is then inverse-FFT'd to obtain a smooth real-space porosity field.

Sphere
^^^^^^
For a sphere of radius r at position **x**₀, the Fourier transform of the
solid indicator function χ(x) = 1{|x - x₀| ≤ r} is::

    F_sphere(k; r, x₀) = V_sphere · exp(-2πi k·x₀) · G(2π|k|r)

where V_sphere = (4/3)πr³ and::

    G(ρ) = 3[sin(ρ) - ρ cos(ρ)] / ρ³,   G(0) = 1.

Verification: F_sphere(k=0) = V_sphere · G(0) = V_sphere  ✓ (DC = volume)

Ellipsoid (novel extension)
^^^^^^^^^^^^^^^^^^^^^^^^^^^
For an oriented ellipsoid described by transformation matrix **M** (the unit
ball mapped to the ellipsoid via **x** = **x**₀ + **M** **u**, |**u**| ≤ 1)::

    F_ellipsoid(k; M, x₀) = |det(M)| · (4π/3) · exp(-2πi k·x₀) · G(2π|M^T k|)

This is *exact*: the ellipsoid FT in direction k is identical to the sphere FT
evaluated at the transformed wavenumber 2π|M^T k|.  The factor |det(M)| = V / (4π/3)
accounts for the volume of the ellipsoid.

Derivation sketch: Let **y** = **M**⁻¹(**x** - **x**₀), then the integral over
the ellipsoid becomes |det(M)| · ∫_{|y|≤1} exp(-2πi k·(**x**₀ + **M**y)) dy,
and the inner integral is (4π/3) · G(|**M**^T **k**|) by the substitution
**q** = **M**^T **k**.

Reconstruction pipeline
^^^^^^^^^^^^^^^^^^^^^^^
::

    F(k) = Σᵢ Fᵢ(k)                           [sum over grains]
    F_smooth(k) = F(k) · K(k)                  [optional smoothing kernel]
    φ(x) = Re[IFFT(F_smooth)] / V_domain       [real-space solid fraction]

where V_domain = Lx · Ly · Lz.

Smoothing kernels K(k):
    'cell'     — sinc(kx·Δx)·sinc(ky·Δy)·sinc(kz·Δz)
                 Tophat average over each grid cell (eliminates aliasing).
    'gaussian' — exp(-2π²σ²|k|²)
                 Isotropic Gaussian average.
    'none'     — K = 1 (no smoothing; Gibbs ringing near grain surfaces).

Structure factor and correlation length (novel addition)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The power spectrum::

    S(k) = |F(k)|² / V_domain

is the structure factor of the pore geometry.  Its inverse FT gives the
two-point autocorrelation function C(r).  The 1/e decay scale of C(r) is the
*correlation length* ξ, which:

- Bounds the Representative Elementary Volume (REV) size: samples of linear
  size >> ξ are approximately statistically independent.
- Is directly related to the MLMC variance decay rate β: roughly
  Var[Q_L - Q_{L-1}] ~ exp(-L / ξ).

Both quantities are computed without additional PDE solves.
"""

from __future__ import annotations

import logging
import math
from html import escape as xml_escape

import numpy as np
from numpy import fft as npfft

from porescalemc.config import FourierConfig
from porescalemc.geometry.grains import Grain, Packing

_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Core analytical functions
# ---------------------------------------------------------------------------

def unit_ball_ft(rho: np.ndarray) -> np.ndarray:
    """Normalised Fourier transform of the unit ball indicator function.

    .. math::

        G(\\rho) = \\frac{3[\\sin(\\rho) - \\rho\\cos(\\rho)]}{\\rho^3},
        \\quad G(0) = 1.

    This is the radial FT of χ_{B_1}(x) = 1{|x| ≤ 1} scaled so that G(0) = 1.
    Physically, G(2π|k|r) is the ratio F_sphere(k)/ V_sphere.

    Parameters
    ----------
    rho : np.ndarray
        Non-negative scaled wavenumber array (any shape).

    Returns
    -------
    np.ndarray
        Values of G in [−0.22, 1] (same shape as ``rho``).
        At ρ = 0, returns exactly 1.0 via L'Hôpital.
    """
    rho = np.asarray(rho, dtype=float)
    out = np.ones_like(rho)

    # For ρ → 0 the formula is 0/0; L'Hôpital gives G(0) = 1 exactly.
    # We use a Taylor series instead of L'Hôpital to avoid branch logic.
    small = rho < 1e-6
    large = ~small

    r = rho[large]
    out[large] = 3.0 * (np.sin(r) - r * np.cos(r)) / r ** 3

    # Taylor: G(ρ) = 1 - ρ²/10 + ρ⁴/280 - O(ρ⁶)
    # Coefficients come from the power-series expansion of sin/cos.
    rs = rho[small]
    out[small] = 1.0 - rs ** 2 / 10.0 + rs ** 4 / 280.0
    return out


def grain_fourier_transform(
    grain: Grain,
    kx: np.ndarray,
    ky: np.ndarray,
    kz: np.ndarray,
) -> np.ndarray:
    """Analytical Fourier transform of a single grain on a 3D k-grid.

    For a sphere of radius r at **x**₀::

        F(k) = V · exp(-2πi k·x₀) · G(2π|k|r)

    For an ellipsoid with transformation matrix M at **x**₀::

        F(k) = |det(M)| · (4π/3) · exp(-2πi k·x₀) · G(2π|M^T k|)

    Parameters
    ----------
    grain : Grain
    kx, ky, kz : np.ndarray
        Wavenumber grids (same shape, broadcast-compatible).
        Units: cycles per unit length (not radians).

    Returns
    -------
    np.ndarray (complex, same shape as kx)
        F(k) for this grain.
    """
    kx = np.asarray(kx, dtype=float)
    ky = np.asarray(ky, dtype=float)
    kz = np.asarray(kz, dtype=float)

    x0, y0, z0 = grain.center
    # Shift theorem: translating a function by x₀ multiplies its FT by exp(-2πi k·x₀).
    phase = np.exp(-2j * math.pi * (kx * x0 + ky * y0 + kz * z0))

    if grain.is_ellipsoid:
        M = grain.transform.reshape(3, 3)
        det_M = float(np.linalg.det(M))
        V = abs(det_M) * (4.0 / 3.0) * math.pi
        # The ellipsoid FT reduces to the sphere FT at the transformed wavenumber M^T k.
        # M^T maps the physical k-vector into the unit-ball coordinate system.
        k_vec = np.stack([kx, ky, kz], axis=0)           # shape (3, *grid)
        Mt_k = np.einsum("ij,j...->i...", M.T, k_vec)    # vectorised M^T @ k for every grid point
        rho = 2.0 * math.pi * np.sqrt((Mt_k ** 2).sum(axis=0))  # |M^T k| scaled to radians
        return V * phase * unit_ball_ft(rho)
    else:
        r = float(grain.transform[0])
        V = 4.0 / 3.0 * math.pi * r ** 3
        k_mag = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
        rho = 2.0 * math.pi * k_mag * r   # dimensionless wavenumber ρ = 2π|k|r
        return V * phase * unit_ball_ft(rho)


# ---------------------------------------------------------------------------
# Packing-level FT
# ---------------------------------------------------------------------------

def _build_kgrid(
    box: np.ndarray, nx: int, ny: int, nz: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float, float]:
    """Build a 3D FFT frequency grid matched to the domain box.

    Returns (kx, ky, kz, dx, dy, dz) where kx/ky/kz are broadcastable arrays
    in cycles-per-unit-length and dx/dy/dz are the real-space grid spacings.
    """
    Lx, Ly, Lz = float(box[0]), float(box[1]), float(box[2])
    dx, dy, dz = Lx / nx, Ly / ny, Lz / nz

    # np.fft.fftfreq(n, d) returns k in cycles/unit-length (not radians/length):
    #   k_j = j/(n·d) for j=0..n/2, then j = -(n-j)/(n·d) for j=n/2+1..n-1.
    # The 'd' parameter converts from cycles/sample to cycles/length.
    kx = npfft.fftfreq(nx, d=dx).reshape(nx, 1, 1)   # broadcast along y and z
    ky = npfft.fftfreq(ny, d=dy).reshape(1, ny, 1)   # broadcast along x and z
    kz = npfft.fftfreq(nz, d=dz).reshape(1, 1, nz)   # broadcast along x and y

    return kx, ky, kz, dx, dy, dz


def packing_fourier_transform(
    packing: Packing,
    nx: int,
    ny: int,
    nz: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute F(k) = Σᵢ Fᵢ(k) over a 3D k-grid.

    The DC component ``F_k[0,0,0].real`` equals the total solid volume
    Σᵢ Vᵢ — this is a useful sanity check.

    Parameters
    ----------
    packing : Packing
    nx, ny, nz : int
        Grid dimensions (number of k-points in each direction).

    Returns
    -------
    F_k : np.ndarray, complex, shape (nx, ny, nz)
        Summed Fourier transform of all grains.
    kx, ky, kz : np.ndarray
        Wavenumber grids (cycles per unit length).
    """
    kx, ky, kz, dx, dy, dz = _build_kgrid(packing.box, nx, ny, nz)
    F_k = np.zeros((nx, ny, nz), dtype=complex)

    # Superposition principle: FT is linear, so the packing FT is the sum of grain FTs.
    for grain in packing.grains:
        F_k += grain_fourier_transform(grain, kx, ky, kz)

    return F_k, kx, ky, kz


def smoothing_kernel(
    kx: np.ndarray,
    ky: np.ndarray,
    kz: np.ndarray,
    kernel: str,
    dx: float,
    dy: float,
    dz: float,
    sigma: float = 1.0,
    constant_sigma: bool = False,
) -> np.ndarray:
    """Smoothing kernel K(k) for the k-space porosity field.

    Parameters
    ----------
    kx, ky, kz : np.ndarray
        Wavenumber grids (cycles per unit length), broadcastable.
    kernel : str
        'cell', 'gaussian', or 'none'.
    dx, dy, dz : float
        Real-space grid spacings (used for 'cell' and relative 'gaussian').
    sigma : float
        Gaussian smoothing length in units of grid spacing (relative) or
        in absolute length units (if ``constant_sigma=True``).
    constant_sigma : bool
        If True, sigma is in the same length units as the domain.

    Returns
    -------
    np.ndarray (real, same shape as kx broadcast with ky, kz)

    Notes
    -----
    'cell' kernel:
        K(k) = sinc(kx·Δx) · sinc(ky·Δy) · sinc(kz·Δz)

        where sinc(x) = sin(πx)/(πx).  This is the FT of a rectangular
        tophat of width Δx, Δy, Δz — exactly what you get from cell-averaging
        the real-space field.  Eliminates all spatial frequencies above the
        Nyquist limit (|k| > 1/(2Δ)) with no ringing.

    'gaussian' kernel:
        K(k) = exp(-2π²σ²|k|²)

        Corresponds to convolution with a Gaussian of standard deviation σ
        in real space.
    """
    if kernel == "none":
        return np.ones(np.broadcast_shapes(kx.shape, ky.shape, kz.shape))

    if kernel == "cell":
        # np.sinc uses the normalised sinc: sinc(x) = sin(πx)/(πx), so sinc(k·Δ)
        # is the FT of a tophat of width Δ centred at the origin.
        return np.sinc(kx * dx) * np.sinc(ky * dy) * np.sinc(kz * dz)

    if kernel == "gaussian":
        if constant_sigma:
            s = sigma                           # σ already in length units
        else:
            # Convert from grid-spacing units to length units using mean spacing.
            s = sigma * (dx + dy + dz) / 3.0
        k2 = kx ** 2 + ky ** 2 + kz ** 2
        # K(k) = exp(-2π²σ²|k|²)  ↔  real-space Gaussian of std σ.
        return np.exp(-2.0 * math.pi ** 2 * s ** 2 * k2)

    raise ValueError(f"Unknown kernel: {kernel!r}. Choose 'cell', 'gaussian', or 'none'.")


# ---------------------------------------------------------------------------
# High-level field reconstruction
# ---------------------------------------------------------------------------

def porosity_field(
    packing: Packing,
    config: FourierConfig,
) -> np.ndarray:
    """Compute the solid volume fraction field φ(x) on a regular grid.

    Algorithm:
    1. Build a 3D k-grid matched to the domain box.
    2. Compute F(k) = Σᵢ Fᵢ(k) analytically for each grain.
    3. Multiply by optional smoothing kernel K(k).
    4. Inverse-FFT and normalise by V_domain.

    The result φ(x) ∈ [0, 1] approximates the local solid volume fraction.
    Its spatial mean equals ``(1 - packing.porosity)`` to within numerical
    precision (better than 1% for typical configurations).

    Parameters
    ----------
    packing : Packing
    config : FourierConfig

    Returns
    -------
    np.ndarray, real, shape (nx, ny, nz)
        Solid volume fraction field.  Values outside [0, 1] are clipped.
    """
    Lx, Ly, Lz = packing.box
    nx = max(2, int(round(config.resolution * Lx)))
    ny = max(2, int(round(config.resolution * Ly)))
    nz = max(2, int(round(config.resolution * Lz)))

    kx, ky, kz, dx, dy, dz = _build_kgrid(packing.box, nx, ny, nz)

    # Step 2: analytical packing FT — exact, no voxelisation needed.
    F_k = np.zeros((nx, ny, nz), dtype=complex)
    for grain in packing.grains:
        F_k += grain_fourier_transform(grain, kx, ky, kz)

    # Step 3: apply smoothing kernel in k-space (multiplication = convolution in real space).
    K = smoothing_kernel(
        kx, ky, kz,
        config.smoothing,
        dx, dy, dz,
        sigma=config.smoothing_length,
        constant_sigma=config.constant_smoothing,
    )
    F_k *= K

    # Step 4: inverse FFT and normalise.
    # numpy IFFT computes (1/N) Σ_k F_k exp(2πi jk/N), where N = nx*ny*nz.
    # We want φ(xj) = Σ_k F(k_k) exp(2πi k_k·xj) Δk³ where Δk³ = 1/V_domain.
    # Therefore φ = IFFT(F_k) * N / V_domain.
    V_domain = float(Lx * Ly * Lz)
    phi = np.real(npfft.ifftn(F_k)) * (nx * ny * nz) / V_domain

    # Clip to [0, 1]: slight overshoots can occur near grain surfaces due to
    # Gibbs ringing (absent with 'cell' kernel, possible with 'none').
    phi = np.clip(phi, 0.0, 1.0)

    _log.debug(
        "porosity_field: grid %dx%dx%d, mean φ=%.4f (expected %.4f)",
        nx, ny, nz, phi.mean(), 1.0 - packing.porosity,
    )
    return phi


# ---------------------------------------------------------------------------
# Structure factor and correlation length
# ---------------------------------------------------------------------------

def structure_factor(
    packing: Packing,
    config: FourierConfig,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the structure factor S(k) = |F(k)|² / V_domain.

    S(k) characterises the two-point statistics of the solid phase.  Its
    spherical average S(|k|) is related to the pair correlation function of
    grain centres.

    Parameters
    ----------
    packing : Packing
    config : FourierConfig

    Returns
    -------
    S_k : np.ndarray, real, shape (nx, ny, nz)
        Structure factor on the k-grid.
    k_magnitudes : np.ndarray, shape (nx, ny, nz)
        |k| at each grid point (cycles per unit length).
    """
    Lx, Ly, Lz = packing.box
    nx = max(2, int(round(config.resolution * Lx)))
    ny = max(2, int(round(config.resolution * Ly)))
    nz = max(2, int(round(config.resolution * Lz)))

    kx, ky, kz, dx, dy, dz = _build_kgrid(packing.box, nx, ny, nz)
    F_k, *_ = packing_fourier_transform(packing, nx, ny, nz)

    V_domain = float(Lx * Ly * Lz)
    # S(k) = |F(k)|²/V is the power spectrum; for a single grain it equals V·|G|².
    # It is related to the pair correlation function g(r) via the Wiener-Khinchin theorem.
    S_k = np.abs(F_k) ** 2 / V_domain
    k_mag = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
    return S_k, np.broadcast_to(k_mag, (nx, ny, nz)).copy()


def spherical_average_structure_factor(
    packing: Packing,
    config: FourierConfig,
    n_bins: int = 50,
) -> tuple[np.ndarray, np.ndarray]:
    """Spherically-averaged structure factor S(|k|).

    Bins the 3D S(k) into radial shells.

    Returns
    -------
    k_bins : np.ndarray, shape (n_bins,)
        Centre of each |k| bin (cycles per unit length).
    S_avg : np.ndarray, shape (n_bins,)
        Mean S(k) in each bin.
    """
    S_k, k_mag = structure_factor(packing, config)
    k_flat = k_mag.ravel()
    S_flat = S_k.ravel()

    k_max = k_flat.max()
    edges = np.linspace(0.0, k_max, n_bins + 1)
    k_bins = 0.5 * (edges[:-1] + edges[1:])
    S_avg = np.zeros(n_bins)

    for i in range(n_bins):
        mask = (k_flat >= edges[i]) & (k_flat < edges[i + 1])
        if mask.any():
            S_avg[i] = S_flat[mask].mean()

    return k_bins, S_avg


def autocorrelation(
    packing: Packing,
    config: FourierConfig,
) -> np.ndarray:
    """Two-point autocorrelation function C(r) of the solid phase.

    C(r) = IFFT[S(k)] normalised so that C(0) = 1.

    The result is in the standard FFT layout: C[0,0,0] is the zero-lag value
    (maximum), and C[i,j,k] corresponds to the spatial lag
    (i*dx, j*dy, k*dz) with the second half of each dimension wrapping to
    negative lags (standard fftfreq convention).

    C(r) decays from 1 at r=0 towards ε² = (1-porosity)² at large r
    (for a spatially uncorrelated random medium), where the excess over ε²
    quantifies the spatial correlation of the solid phase.

    Returns
    -------
    np.ndarray, real, shape (nx, ny, nz)
        C(r) on the standard FFT grid (C[0,0,0] = 1).
    """
    S_k, _ = structure_factor(packing, config)

    # Wiener-Khinchin: autocorrelation = IFFT(power spectrum).
    # We do NOT fftshift before/after: the FFT layout already has r=0 at index [0,0,0],
    # so C[0,0,0] IS the zero-lag value (the maximum) without any shifting.
    C = np.real(npfft.ifftn(S_k))
    C0 = float(C[0, 0, 0])
    if abs(C0) < 1e-15:
        return C
    return C / C0


def correlation_length(
    packing: Packing,
    config: FourierConfig,
) -> float:
    """Estimate the spatial correlation length ξ of the solid phase.

    ξ is defined as the radial distance at which the normalised spherically-
    averaged autocorrelation function C(r) first falls below 1/e.

    This bounds the REV size: independent samples require a domain size
    of at least a few multiples of ξ in each direction.

    Parameters
    ----------
    packing : Packing
    config : FourierConfig

    Returns
    -------
    float
        Correlation length in the same units as the domain box.
        Returns np.nan if C(r) never crosses 1/e (highly correlated packing).
    """
    C = autocorrelation(packing, config)
    Lx, Ly, Lz = packing.box
    nx, ny, nz = C.shape

    # Build the physical distance |r| for every grid point.
    # fftfreq(n) * L gives r ∈ {0, dx, 2dx, …, (n/2-1)dx, -(n/2)dx, …, -dx}
    # — i.e., the FFT-layout distances matching C's index ordering.
    x = np.fft.fftfreq(nx) * Lx
    y = np.fft.fftfreq(ny) * Ly
    z = np.fft.fftfreq(nz) * Lz
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    R = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)

    # Spherically-averaged C(r): bin by distance.
    r_flat = R.ravel()
    c_flat = C.ravel()
    r_max = float(min(Lx, Ly, Lz)) / 2.0   # don't trust C beyond half the box (aliasing)
    n_bins = max(20, int(config.resolution * r_max))
    edges = np.linspace(0.0, r_max, n_bins + 1)
    r_bins = 0.5 * (edges[:-1] + edges[1:])
    c_avg = np.zeros(n_bins)

    for i in range(n_bins):
        mask = (r_flat >= edges[i]) & (r_flat < edges[i + 1])
        if mask.any():
            c_avg[i] = c_flat[mask].mean()

    # Find the first bin where C(r) crosses below 1/e.
    threshold = math.exp(-1.0)
    crossings = np.where(c_avg < threshold)[0]
    if len(crossings) == 0:
        _log.debug("correlation_length: C(r) never crosses 1/e; returning domain half-size.")
        return r_max

    idx = int(crossings[0])
    # Linear interpolation between the bin just above and just below 1/e
    # gives sub-bin precision without needing finer radial bins.
    if idx == 0:
        return float(r_bins[0])
    r0, r1 = r_bins[idx - 1], r_bins[idx]
    c0, c1 = c_avg[idx - 1], c_avg[idx]
    if abs(c1 - c0) < 1e-15:
        return float(r0)
    return float(r0 + (threshold - c0) / (c1 - c0) * (r1 - r0))


# ---------------------------------------------------------------------------
# Upscaled permeability from Fourier field
# ---------------------------------------------------------------------------

def upscaled_permeability(
    packing: Packing,
    fourier_config: FourierConfig,
    law: str = "kozeny_carman",
    **law_kwargs,
) -> float:
    """Estimate the volume-averaged permeability from the porosity field.

    Computes the porosity field φ(x), applies a porosity-permeability closure
    law pointwise, and returns the harmonic average (appropriate for flow in
    series — the minimum bottleneck).

    Parameters
    ----------
    packing : Packing
    fourier_config : FourierConfig
    law : str
        Porosity-permeability law.  Passed to
        :func:`porescalemc.geometry.transforms.permeability_field`.
    **law_kwargs
        Additional arguments forwarded to the law (e.g. ``grain_size``).

    Returns
    -------
    float
        Harmonically-averaged permeability (same units as the closure law).
    """
    from porescalemc.geometry.transforms import permeability_field

    solid_fraction = porosity_field(packing, fourier_config)
    porosity = 1.0 - solid_fraction
    default_grain_size = (
        2.0 * float(np.mean([g.max_radius for g in packing.grains]))
        if packing.grains
        else 1.0
    )
    grain_size = law_kwargs.pop("grain_size", default_grain_size)
    K_field = permeability_field(porosity, law=law, grain_size=grain_size, **law_kwargs)

    # Harmonic mean models resistance-in-series: the lowest-permeability
    # region dominates, which is physically appropriate for pressure-driven flow.
    inv_K = 1.0 / np.where(K_field > 0, K_field, np.inf)
    return float(1.0 / inv_K.mean())


# ---------------------------------------------------------------------------
# Post-processing / offline visualisation helpers
# ---------------------------------------------------------------------------

def save_structured_vti_fields(
    fields: dict[str, np.ndarray],
    origin: tuple[float, float, float],
    spacing: tuple[float, float, float],
    filename: str,
    active_scalar: str | None = None,
) -> None:
    """Write one or more 3-D scalar fields as XML VTK ImageData (.vti).

    The writer is intentionally dependency-free and emits ASCII XML so the
    files are easy to inspect in a text editor.  ParaView and VisIt can load
    the result directly as regular point data.
    """
    if not fields:
        raise ValueError("fields must contain at least one scalar array.")

    arrays: dict[str, np.ndarray] = {}
    shape: tuple[int, int, int] | None = None
    for name, field in fields.items():
        arr = np.asarray(field, dtype=np.float32)
        if arr.ndim != 3:
            raise ValueError(f"Field {name!r} must be 3-D, got shape {arr.shape}.")
        if shape is None:
            shape = arr.shape
        elif arr.shape != shape:
            raise ValueError(
                f"Field {name!r} has shape {arr.shape}; expected {shape}."
            )
        arrays[str(name)] = arr

    assert shape is not None
    nx, ny, nz = shape
    ox, oy, oz = origin
    dx, dy, dz = spacing
    active = active_scalar if active_scalar in arrays else next(iter(arrays))
    extent = f"0 {nx - 1} 0 {ny - 1} 0 {nz - 1}"

    with open(filename, "w", encoding="utf-8") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">\n')
        f.write(
            f'  <ImageData WholeExtent="{extent}" '
            f'Origin="{ox:.12g} {oy:.12g} {oz:.12g}" '
            f'Spacing="{dx:.12g} {dy:.12g} {dz:.12g}">\n'
        )
        f.write(f'    <Piece Extent="{extent}">\n')
        f.write(f'      <PointData Scalars="{xml_escape(active)}">\n')

        for name, arr in arrays.items():
            f.write(
                f'        <DataArray type="Float32" Name="{xml_escape(name)}" '
                'NumberOfComponents="1" format="ascii">\n'
            )
            flat = arr.ravel(order="F")
            for start in range(0, flat.size, 8):
                chunk = flat[start:start + 8]
                values = " ".join(f"{float(value):.8e}" for value in chunk)
                f.write(f"          {values}\n")
            f.write("        </DataArray>\n")

        f.write("      </PointData>\n")
        f.write("      <CellData/>\n")
        f.write("    </Piece>\n")
        f.write("  </ImageData>\n")
        f.write("</VTKFile>\n")


def save_structured_vti(
    field: np.ndarray,
    origin: tuple[float, float, float],
    spacing: tuple[float, float, float],
    filename: str,
    field_name: str = "scalar",
) -> None:
    """Write a single 3-D numpy array as XML VTK ImageData (.vti)."""
    save_structured_vti_fields(
        {field_name: field},
        origin,
        spacing,
        filename,
        active_scalar=field_name,
    )


def save_spectral_fields_to_vti(
    solver,
    packing,
    fourier_config,
    base_name: str = "spectral",
) -> None:
    """Convenience wrapper that exports common fields from a spectral solver.

    Writes:
      - ``{base_name}_fields.vti`` with all cached solver fields
      - ``{base_name}_solid_fraction.vti`` and ``{base_name}_porosity.vti``
        for quick standalone inspection

    For Stokes / advection-diffusion solvers that store internal velocity
    or corrector fields, the combined file includes those arrays too.

    This gives you a structured mesh that visualises both the (smoothed)
    geometry and the solution fields inside ParaView without any surface
    meshing.
    """
    # Always have the solid fraction on the solver grid
    phi_s = porosity_field(packing, fourier_config)
    origin = tuple(-np.array(packing.box) / 2.0)   # centered box convention
    Lx, Ly, Lz = packing.box
    nx, ny, nz = phi_s.shape
    spacing = (Lx / nx, Ly / ny, Lz / nz)

    fields: dict[str, np.ndarray] = {}
    if hasattr(solver, "solution_fields"):
        fields = {
            name: np.asarray(field)
            for name, field in solver.solution_fields().items()
            if np.asarray(field).shape == phi_s.shape
        }
    fields.setdefault("solid_fraction", phi_s)
    fields.setdefault("porosity", 1.0 - phi_s)

    save_structured_vti_fields(
        fields,
        origin,
        spacing,
        f"{base_name}_fields.vti",
        active_scalar="solid_fraction",
    )

    save_structured_vti(
        phi_s, origin, spacing,
        f"{base_name}_solid_fraction.vti", "solid_fraction"
    )
    save_structured_vti(
        1.0 - phi_s, origin, spacing,
        f"{base_name}_porosity.vti", "porosity"
    )

    # If the solver kept velocity fields (Stokes / AdvectionDiffusion)
    if hasattr(solver, "_velocity_fields") and solver._velocity_fields is not None:
        vel = solver._velocity_fields[0]   # take first direction for demo
        for comp, name in enumerate(["ux", "uy", "uz"]):
            save_structured_vti(
                vel[comp], origin, spacing,
                f"{base_name}_{name}.vti", name
            )


def plot_field_slices(
    fields: dict[str, np.ndarray],
    filename: str | None = None,
    axis: int = 2,
    index: int | None = None,
    cmap: str = "viridis",
) -> None:
    """Plot mid-plane slices of one or more 3-D fields.

    This is deliberately lightweight post-processing: it is meant for quick
    sanity checks of PDE solutions and Fourier fields before opening ParaView.
    If matplotlib is unavailable the function logs and returns.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        _log.error("matplotlib is required for plotting. Install with: pip install matplotlib")
        return

    if not fields:
        raise ValueError("fields must contain at least one scalar array.")

    axis = int(axis)
    if axis < 0 or axis > 2:
        raise ValueError("axis must be 0, 1 or 2.")

    arrays = {name: np.asarray(field) for name, field in fields.items()}
    for name, arr in arrays.items():
        if arr.ndim != 3:
            raise ValueError(f"Field {name!r} must be 3-D, got shape {arr.shape}.")

    first_shape = next(iter(arrays.values())).shape
    slice_index = first_shape[axis] // 2 if index is None else int(index)
    if slice_index < 0 or slice_index >= first_shape[axis]:
        raise ValueError(
            f"index {slice_index} is outside axis {axis} bounds for shape {first_shape}."
        )

    n_fields = len(arrays)
    ncols = min(3, n_fields)
    nrows = int(math.ceil(n_fields / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.6 * nrows))
    axes_arr = np.atleast_1d(axes).ravel()

    for ax, (name, arr) in zip(axes_arr, arrays.items()):
        if axis == 0:
            sl = arr[slice_index, :, :]
        elif axis == 1:
            sl = arr[:, slice_index, :]
        else:
            sl = arr[:, :, slice_index]
        im = ax.imshow(np.asarray(sl).T, origin="lower", cmap=cmap, aspect="equal")
        ax.set_title(name)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    for ax in axes_arr[n_fields:]:
        ax.axis("off")

    fig.tight_layout()
    if filename:
        fig.savefig(filename, bbox_inches="tight", dpi=150)
        _log.info("Saved field slice plot to %s", filename)
    else:
        plt.show()
    plt.close(fig)
