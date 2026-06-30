"""
Typed configuration dataclasses for porescalemc.

All configuration lives here — no exec(open(...)), no global mutation.
Pass config objects explicitly to every function that needs them.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field


# ---------------------------------------------------------------------------
# Packing / geometry configuration
# ---------------------------------------------------------------------------

@dataclass
class PackingConfig:
    """Parameters that fully describe a random grain packing experiment.

    Grain size distribution (PSD):
        'lognormal' — log-normal with mean ``mu`` and CoV ``coeff_var``
        'uniform'   — uniform on [mu*(1-coeff_var), mu*(1+coeff_var)]
        'constant'  — all grains have radius exactly ``mu``

    Hierarchy flags (used by MLMC level scaling):
        The ``hierarchy`` string in MLMCConfig controls which parameters are
        scaled per level.  This object holds the base (level-0) values.
    """

    # --- Domain geometry ---
    # xlen/ylen/zlen define the rectangular periodic box; dimension is used
    # only by hierarchy 'n' to scale grain count with box volume.
    dimension: int = 3
    xlen: float = 1.0
    ylen: float = 1.0
    zlen: float = 1.0

    # --- Grain size distribution ---
    # mu is the mean radius; coeff_var = std/mean parameterises spread.
    # For lognormal PSD: sigma_ln = sqrt(log(cv²+1)), mu_ln = log(mu) - sigma_ln²/2.
    mu: float = 0.1           # mean grain radius
    n_grains: int = 20        # target number of grains
    coeff_var: float = 0.1    # coefficient of variation of radius distribution
    psd: str = "lognormal"    # 'lognormal' | 'uniform' | 'constant'

    # --- Grain shape ---
    ellipsoid: bool = False   # True → sample Haar-random orientation + 3 semi-axes

    # --- Boundary conditions ---
    # periodic=True: grains wrap across box faces (minimum-image distances).
    # detached=True: reject placements that overlap any already-placed grain.
    # detached_bc: margin multiplier so grains don't touch box walls (0 = disabled).
    periodic: bool = True
    detached: bool = True
    detached_bc: float = 1.0  # grain must stay ≥ detached_bc * r from each face

    # --- RSA stopping criteria ---
    min_porosity: float = 0.3   # abort grain addition below this porosity
    max_tries: int = 1000       # rejection-sampling attempts per grain before giving up

    # --- Extra margin around box ---
    void_space: float = 0.0   # gap added to detached_bc margin (units of mean radius)

    # --- Jodrey-Tory post-processing ---
    # JT iteratively resolves overlaps by moving the closest pairs apart.
    # jt_max_dist < 0 disables the optional attraction phase.
    use_jodrey_tory: bool = False
    jt_min_dist: float = 1.0    # target minimum normalised distance (1 = just touching)
    jt_max_dist: float = -1.0   # attraction threshold (pull far pairs together); ≤0 = off
    jt_n_moves: int = 10        # number of pairs to move per JT iteration
    jt_eps: float = 0.5         # step size fraction of the overlap gap
    jt_cluster: int = 3         # nearest-neighbour count for the attraction phase

    @property
    def box(self):
        """Domain lengths as a tuple (Lx, Ly, Lz)."""
        return (self.xlen, self.ylen, self.zlen)

    @property
    def volume(self) -> float:
        return self.xlen * self.ylen * self.zlen


# ---------------------------------------------------------------------------
# MLMC configuration
# ---------------------------------------------------------------------------

@dataclass
class MLMCConfig:
    """Parameters for the multilevel Monte Carlo estimator.

    Hierarchy flags (``hierarchy`` string):
        'g' — grid refinement: mesh resolution doubles each level
        'd' — domain scaling: domain grows each level
        's' — stone (grain) size: mean radius shrinks each level
        'n' — number of grains: grain count increases each level

    Multiple flags can be combined, e.g. ``hierarchy='gn'``.

    MLMC pair format (see mlmc/statistics.py):
        At each level l, each sample is the concatenated vector
        ``[Q_l - Q_{l-1}, Q_l]``.  For l=0: ``[Q_0, Q_0]``.
    """

    # --- Level structure ---
    n_levels: int = 3          # initial number of levels
    max_levels: int = 5        # adaptive algorithm cannot exceed this

    # --- Sample counts ---
    # mratio controls how many FEWER samples we draw at each finer level:
    #   M_l ≈ M_0 / mratio^l  (more expensive levels get fewer samples).
    # This is independent of refratio, which governs the geometry resolution.
    m0: int = 32               # samples at coarsest level (level 0)
    mratio: float = 4.0        # sample-count ratio between adjacent levels
    min_samples: int = 5       # floor per level for numerical stability of variance

    # --- Spatial refinement ---
    # refratio controls how much finer each level's geometry/mesh is.
    # Classic MLMC uses refratio=2 (halve grid spacing each level).
    # Kept separate from mratio because optimal M_l depends on Var_l/W_l,
    # not on the mesh refinement factor.
    refratio: float = 2.0      # spatial refinement ratio (dx_l = dx_0 / refratio^l)

    # --- Convergence criteria ---
    tolerance: float = 1e-2   # overall relative RMSE target TOL
    error_split: float = 0.5   # θ: fraction of TOL assigned to bias (1-θ → stat error)
    confidence: float = 0.99   # confidence level for the CLT-based stat error bound

    # --- Execution ---
    n_workers: int = 1         # parallel processes (1 = serial, avoids pickle overhead)
    algorithm: str = "fixed"   # 'fixed' | 'adaptive' (Giles 2015 adaptive algorithm)
    estimator_type: str = "pair"  # 'pair' (standard) | 'triplet' (future)

    # --- Hierarchy type ---
    hierarchy: str = "g"       # see class docstring

    reuse_samples: bool = True  # accumulate samples across adaptive iterations

    # --- Assumed convergence rates ---
    # These three form the MLMC complexity triangle:
    #   α (weak)   : |E[Q_l] - Q*| ~ C_α · refratio^{-α·l}  (bias decay)
    #   β (strong) : Var[Q_l - Q_{l-1}] ~ C_β · refratio^{-β·l}  (variance decay)
    #   γ (work)   : W_l ~ C_γ · refratio^{γ·l}  (cost growth)
    # Complexity theorem: total work ~ O(TOL^{-2}) when β > γ/α.
    beta: float = 1.2    # variance decay rate
    alpha: float = 2.5   # bias decay rate
    gamma: float = 3.0   # work growth rate

    # --- Bias correction ---
    extrapolate_bias: bool = False  # add Richardson correction: +means[-1]/(2^α-1)
    bias_mode: str = "last_increment"  # 'last_increment' | 'adaptive'

    @property
    def calpha(self) -> float:
        """Normal quantile z_{alpha} for the confidence level."""
        from scipy.stats import norm
        return float(norm.ppf(self.confidence))

    @property
    def initial_sample_counts(self) -> list[int]:
        """Initial number of samples per level [M_0, M_1, ..., M_{L-1}]."""
        counts = [int(self.m0)]
        for _ in range(self.n_levels - 1):
            counts.append(max(self.min_samples, int(counts[-1] / self.mratio)))
        return counts


# ---------------------------------------------------------------------------
# Fourier field configuration
# ---------------------------------------------------------------------------

@dataclass
class FourierConfig:
    """Parameters for the analytical Fourier-space porosity field.

    The porosity field is reconstructed from the analytical FT of each grain
    (see geometry/fourier_field.py for the full mathematical derivation).

    Smoothing kernels:
        'cell'     — sinc(kx·Δx)·sinc(ky·Δy)·sinc(kz·Δz)
                     Corresponds to a tophat average over each grid cell.
                     Eliminates all spatial frequencies above the Nyquist limit.
        'gaussian' — exp(−2π²σ²|k|²)
                     Isotropic Gaussian average.  ``smoothing_length`` sets σ.
        'none'     — no smoothing; Gibbs ringing may appear near grain boundaries.
    """

    resolution: float = 10.0      # grid points per unit length per dimension
    smoothing: str = "cell"        # 'cell' | 'gaussian' | 'none'
    smoothing_length: float = 1.0  # σ in units of grid spacing (gaussian only)
    constant_smoothing: bool = False  # if True, σ is in absolute length units


# ---------------------------------------------------------------------------
# Hierarchy-level helper
# ---------------------------------------------------------------------------

@dataclass
class HierarchyLevel:
    """Pre-computed scaling factors for a single MLMC level."""

    level: int
    mesh_refinement_factor: float  # for 'g': dx = dx0 / refratio^level
    domain_scale: float            # for 'd': Lx = Lx0 * domain_scale
    grain_size_scale: float        # for 's': mu = mu0 * grain_size_scale
    n_grains_scale: float          # for 'n': n = n0 * n_grains_scale


def hierarchy_level(
    mlmc_cfg: MLMCConfig,
    packing_cfg: PackingConfig,
    level: int,
) -> tuple[PackingConfig, float]:
    """Return a scaled PackingConfig and mesh-resolution factor for ``level``.

    Scaling rules (per flag in mlmc_cfg.hierarchy):
        'g' — mesh resolution multiplied by refratio^level (PackingConfig unchanged)
        'd' — domain lengths scaled by refratio^level (larger domain = lower bias for 'd')
        's' — mean grain radius divided by refratio^level (finer grains = finer scale)
        'n' — grain count multiplied by refratio^(dimension*level) (constant packing fraction)

    Returns
    -------
    packing : PackingConfig
        Scaled config for this level.
    mesh_factor : float
        Mesh resolution multiplier for this level (refratio^level for 'g', else 1).
    """
    import copy
    p = copy.deepcopy(packing_cfg)
    h = mlmc_cfg.hierarchy
    # Use refratio (spatial), NOT mratio (sample count) — these are independent.
    ref = mlmc_cfg.refratio

    if "d" in h:
        # Domain grows by refratio each level; keeps grain density constant.
        scale = ref ** level
        p.xlen *= scale
        p.ylen *= scale
        p.zlen *= scale

    if "s" in h:
        # Finer grains at higher levels → better resolved geometry.
        p.mu /= ref ** level

    if "n" in h:
        # Scale grain count so solid fraction stays constant:
        # same density in a larger domain means more grains.
        p.n_grains = int(packing_cfg.n_grains * (ref ** level) ** packing_cfg.dimension)

    # mesh_factor is passed to Fourier/mesh solvers to set their grid resolution.
    mesh_factor = ref ** level if "g" in h else 1.0

    return p, mesh_factor
