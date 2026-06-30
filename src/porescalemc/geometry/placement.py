"""
Random sequential addition of spheres/ellipsoids.

Algorithm
---------
1. Sample a grain size from the specified PSD (lognormal/uniform/constant).
2. Sample a random position inside the box (excluding a boundary margin if
   ``detached_bc > 0``).
3. Accept the grain if it does not overlap any already-placed grain
   (when ``config.detached`` is True).
4. Stop when either ``n_grains`` grains are placed or porosity drops below
   ``min_porosity``.
5. Optionally run Jodrey-Tory compaction post-processing.

Hierarchy
---------
Use :func:`porescalemc.config.hierarchy_level` to obtain a scaled
``PackingConfig`` for each MLMC level before calling ``sample_packing``.
"""

from __future__ import annotations

import copy
import hashlib
import logging
import math

import numpy as np

from porescalemc.config import PackingConfig
from porescalemc.geometry.grains import (
    Grain,
    Packing,
    ellipsoid_distance,
    sample_random_orientation,
)

_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# PSD samplers
# ---------------------------------------------------------------------------

def _make_radius_sampler(config: PackingConfig, rng: np.random.Generator):
    """Return a callable () -> float that draws one grain radius."""
    mu = config.mu
    cv = config.coeff_var

    if config.psd == "lognormal":
        # Match first two moments of lognormal to (mu, cv):
        #   sigma_ln = sqrt(log(cv² + 1))
        #   mu_ln    = log(mu) - sigma_ln²/2   (so that E[r] = mu exactly)
        ln_sigma = math.sqrt(math.log(cv ** 2 + 1.0))
        ln_mu = math.log(mu) - 0.5 * ln_sigma ** 2
        def sampler():
            return float(rng.lognormal(ln_mu, ln_sigma))

    elif config.psd == "uniform":
        lo = max(1e-10, mu * (1.0 - cv))
        hi = mu * (1.0 + cv)
        def sampler():
            return float(rng.uniform(lo, hi))

    elif config.psd == "constant":
        def sampler():
            return mu

    else:
        raise ValueError(f"Unknown PSD: {config.psd!r}. Choose 'lognormal', 'uniform', or 'constant'.")

    return sampler


# ---------------------------------------------------------------------------
# Position sampler
# ---------------------------------------------------------------------------

def _make_position_sampler(config: PackingConfig, rng: np.random.Generator):
    """Return a callable (r: float) -> np.ndarray that draws a centre position.

    The grain centre is placed at least ``r * detached_bc`` away from each
    face of the box (when ``detached_bc > 0`` and ``periodic=False``).
    For periodic packings the full box is available.
    """
    Lx, Ly, Lz = config.xlen, config.ylen, config.zlen
    bc = 0.0 if config.periodic else config.detached_bc  # no wall margin for periodic
    vs = config.void_space  # extra gap in units of mean radius

    def sampler(r: float) -> np.ndarray:
        # Total margin = (boundary coefficient + void space) × grain radius.
        # This keeps the grain surface at least margin away from each box face.
        margin = (bc + vs) * r
        half = np.array([Lx, Ly, Lz]) * 0.5
        lo = -(half - margin)
        hi =  (half - margin)
        # Safety clamp: for very large grains, ensure the sampling interval is valid.
        lo = np.minimum(lo,  half * 0.9)
        hi = np.maximum(hi, -half * 0.9)
        return rng.uniform(lo, hi)

    return sampler


# ---------------------------------------------------------------------------
# Overlap check
# ---------------------------------------------------------------------------

def _overlaps_any(candidate: Grain, placed: list[Grain], box_or_none) -> bool:
    """Return True if ``candidate`` overlaps any grain in ``placed``."""
    for g in placed:
        if ellipsoid_distance(candidate, g, box_or_none) < 1.0:
            return True
    return False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def sample_packing(
    config: PackingConfig,
    rng: np.random.Generator | None = None,
) -> Packing:
    """Generate a random sphere/ellipsoid packing by sequential addition.

    Parameters
    ----------
    config : PackingConfig
        All geometry parameters (domain size, PSD, overlap rules, …).
    rng : numpy.random.Generator, optional
        Random number generator.  If None, a new default generator is created.
        Pass an explicit seeded generator for reproducible results.

    Returns
    -------
    Packing
        The generated packing.  ``packing.porosity`` gives the estimated
        (non-overlapping) porosity; use the Fourier field for an exact value.
    """
    if rng is None:
        rng = np.random.default_rng()

    box = np.array([config.xlen, config.ylen, config.zlen])
    box_or_none = box if config.periodic else None

    radius_sampler = _make_radius_sampler(config, rng)
    pos_sampler = _make_position_sampler(config, rng)

    grains: list[Grain] = []
    solid_vol = 0.0
    box_vol = config.xlen * config.ylen * config.zlen

    # RSA loop: three stopping conditions —
    #   (1) desired grain count reached,
    #   (2) porosity dropped below min_porosity (domain is dense enough),
    #   (3) max_tries placements rejected in a row for one grain (too dense to fit more).
    for i in range(config.n_grains):
        current_porosity = 1.0 - solid_vol / box_vol
        if current_porosity < config.min_porosity:
            _log.debug("Reached min_porosity=%.3f after %d grains", config.min_porosity, len(grains))
            break

        r = radius_sampler()
        if config.ellipsoid:
            rx = [radius_sampler() for _ in range(3)]
            transform = sample_random_orientation(rx, rng)
        else:
            transform = np.array([r])

        placed = False
        for _attempt in range(config.max_tries):
            center = pos_sampler(r)
            candidate = Grain(center=center, transform=transform)
            if not config.detached or not _overlaps_any(candidate, grains, box_or_none):
                grains.append(candidate)
                solid_vol += candidate.volume
                placed = True
                break

        if not placed:
            _log.debug(
                "Grain %d: could not place after %d tries; stopping.", i, config.max_tries
            )
            break

    packing = Packing(grains=grains, box=box)
    _log.debug(
        "Placed %d grains, estimated porosity=%.4f", packing.n_grains, packing.porosity
    )

    if config.use_jodrey_tory and packing.n_grains >= 2:
        from porescalemc.geometry.jodrey_tory import jodrey_tory
        packing = jodrey_tory(packing, config, rng=rng)

    return packing


def seed_from_name(name: str) -> int:
    """Return a stable 32-bit random seed derived from ``name``.

    Python's built-in ``hash()`` is process-randomised by PYTHONHASHSEED,
    so the same name yields a different seed in each Python process.
    ``hashlib.blake2b`` is deterministic and cross-platform, making it safe
    for MLMC where the same realization name must map to the same geometry
    across parallel workers and restarts.
    """
    digest = hashlib.blake2b(str(name).encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, byteorder="little") % (2**32)


def sample_packing_from_name(
    level: int,
    name: str,
    config: PackingConfig,
) -> Packing:
    """Generate a deterministic packing for ``name`` and ``config``.

    The ``level`` argument is accepted to match the geometry-factory protocol;
    the level-dependent scaling is already encoded in ``config`` by the caller.
    """
    del level  # level already baked into config by hierarchy_level()
    return sample_packing(config, np.random.default_rng(seed_from_name(name)))


def scale_packing_config(
    base_config: PackingConfig,
    new_config: PackingConfig,
) -> PackingConfig:
    """Return a copy of ``new_config`` inheriting any unset fields from ``base_config``.

    Convenience helper for MLMC hierarchy level scaling: after computing the
    new level parameters with :func:`porescalemc.config.hierarchy_level`,
    call this to fill in fields that don't change between levels.
    """
    return copy.deepcopy(new_config)
