"""
Jodrey-Tory (JT) packing compaction algorithm.

The JT algorithm iteratively resolves grain–grain overlaps by moving each
overlapping pair apart along their centre-to-centre vector.  It can also
optionally push grains that are too far apart closer together (controlled by
``config.jt_max_dist``).

Reference: Jodrey & Tory, Physical Review A 32 (1985) 2347.

This implementation:
- Works for both spheres and oriented ellipsoids (using the same normalised
  distance metric as :func:`porescalemc.geometry.grains.ellipsoid_distance`).
- Supports periodic boundary conditions via minimum-image convention.
- Returns a new Packing; the input is not mutated.
"""

from __future__ import annotations

import logging

import numpy as np

from porescalemc.config import PackingConfig
from porescalemc.geometry.grains import (
    Grain,
    Packing,
    _periodic_displacement,
    ellipsoid_distance,
    stretching_length,
)

_log = logging.getLogger(__name__)


def jodrey_tory(
    packing: Packing,
    config: PackingConfig,
    rng: np.random.Generator | None = None,
) -> Packing:
    """Run the Jodrey-Tory compaction algorithm on ``packing``.

    Iteratively moves overlapping grain pairs apart until the minimum
    normalised distance between any pair exceeds ``config.jt_min_dist``,
    or ``config.max_tries`` iterations are reached.

    The algorithm has two phases:
    - **Repulsion phase** (always active): find the n_moves most-overlapping
      pairs (d < jt_min_dist) and push each pair apart by eps * overlap_gap.
    - **Attraction phase** (optional, when jt_max_dist > 0): once all overlaps
      are resolved, pull the farthest pairs among each grain's neighbourhood
      together to achieve a target compaction.

    Parameters
    ----------
    packing : Packing
        Input packing.  Not mutated.
    config : PackingConfig
        Uses: ``jt_min_dist``, ``jt_max_dist``, ``jt_n_moves``, ``jt_eps``,
        ``jt_cluster``, ``max_tries``, ``periodic``.
    rng : numpy.random.Generator, optional
        Unused currently; reserved for future stochastic perturbations.

    Returns
    -------
    Packing
        New packing with reduced overlaps.
    """
    n = packing.n_grains
    if n < 2:
        return packing

    centers = packing.centers.copy()  # (N, 3) — mutable working copy
    grains = packing.grains
    box = packing.box
    box_or_none = box if config.periodic else None

    min_dist = config.jt_min_dist
    max_dist = config.jt_max_dist
    eps = config.jt_eps
    n_moves = config.jt_n_moves
    cluster = config.jt_cluster

    def _rebuild_grains():
        """Reconstruct Grain list from the current mutable centers array."""
        return [Grain(center=centers[i].copy(), transform=grains[i].transform) for i in range(n)]

    def _dist_matrix():
        """Compute pairwise normalised distances and their grain indices."""
        pairs = []
        pair_indices = []
        for i in range(n):
            for j in range(i + 1, n):
                g_i = Grain(center=centers[i], transform=grains[i].transform)
                g_j = Grain(center=centers[j], transform=grains[j].transform)
                pairs.append(ellipsoid_distance(g_i, g_j, box_or_none))
                pair_indices.append((i, j))
        return np.array(pairs), pair_indices

    def _apply_boundary(indices_and_grains):
        """Clip or wrap moved grains back into the simulation box."""
        if config.periodic:
            for idx, _grain in indices_and_grains:
                for dim in range(3):
                    L = box[dim]
                    # Wrap: shift by L/2, apply modulo, then shift back.
                    centers[idx, dim] = ((centers[idx, dim] + L / 2) % L) - L / 2
        else:
            bc = config.detached_bc
            for idx, grain in indices_and_grains:
                for dim in range(3):
                    r_dim = stretching_length(grain, np.eye(3)[dim])
                    limit = max(0.0, box[dim] / 2.0 - bc * r_dim)
                    centers[idx, dim] = float(np.clip(centers[idx, dim], -limit, limit))

    def _move_pair(i: int, j: int, d_norm: float, target: float, sign: float) -> None:
        """Move a pair toward (sign=-1) or away from (sign=+1) each other.

        Each grain moves by half the displacement, keeping the centre of mass fixed.
        Displacement = eps * |target - d_norm| * denom, where denom is the sum of
        stretching lengths (the physical "diameter" in the contact direction).
        """
        g_i = Grain(center=centers[i], transform=grains[i].transform)
        g_j = Grain(center=centers[j], transform=grains[j].transform)

        if box_or_none is not None:
            dvec = _periodic_displacement(centers[i], centers[j], box)
        else:
            dvec = centers[j] - centers[i]

        dist_eucl = float(np.linalg.norm(dvec))
        if dist_eucl < 1e-12:
            # Coincident grains: perturb in a random direction to break degeneracy.
            dvec = rng.standard_normal(3) if rng is not None else np.ones(3)
            dist_eucl = float(np.linalg.norm(dvec))

        direction = dvec / dist_eucl
        # denom = sum of stretching lengths in the contact direction.
        denom = stretching_length(g_i, direction) + stretching_length(g_j, direction)
        if denom < 1e-15:
            return

        # Step size is proportional to how far we are from the target distance.
        delta = eps * abs(target - d_norm) * denom
        # Move both grains symmetrically; sign=+1 pushes apart, sign=-1 pulls together.
        centers[i] -= sign * 0.5 * delta * direction
        centers[j] += sign * 0.5 * delta * direction
        _apply_boundary([(i, g_i), (j, g_j)])

    for iteration in range(config.max_tries):
        dist_flat, pair_indices = _dist_matrix()
        if len(dist_flat) == 0:
            break

        # --- Repulsion phase: resolve overlaps (d < jt_min_dist). ---
        overlapping = [
            (idx, pair_indices[idx], dist)
            for idx, dist in enumerate(dist_flat)
            if dist < min_dist
        ]

        if overlapping:
            # Move the n_moves most-overlapping pairs (smallest d first).
            selected = sorted(overlapping, key=lambda item: item[2])[:n_moves]
            for _idx, (i, j), d_norm in selected:
                _move_pair(i, j, d_norm, min_dist, sign=1.0)
            continue  # skip attraction phase while overlaps remain

        if max_dist <= 0:
            # No attraction phase requested; we're done as soon as overlaps are gone.
            _log.debug("JT converged at iteration %d (min_dist=%.4f)", iteration, dist_flat.min())
            break

        # --- Attraction phase: pull far neighbours together (d > jt_max_dist). ---
        # For each grain, look at its `cluster` nearest neighbours and attract any
        # that are farther than jt_max_dist.  This raises solid fraction without
        # re-introducing overlaps.
        candidates = {}
        for i in range(n):
            local = []
            for idx, (a, b) in enumerate(pair_indices):
                if a == i or b == i:
                    local.append((dist_flat[idx], idx, (a, b)))
            local.sort(key=lambda item: item[0])
            for d_norm, idx, pair in local[:cluster]:
                if d_norm > max_dist:
                    candidates[idx] = (idx, pair, d_norm)

        if not candidates:
            _log.debug("JT converged at iteration %d with max-distance phase.", iteration)
            break

        selected = sorted(candidates.values(), key=lambda item: item[2], reverse=True)[:n_moves]
        for _idx, (i, j), d_norm in selected:
            _move_pair(i, j, d_norm, max_dist, sign=-1.0)
    else:
        _log.warning("JT: reached max_tries=%d without converging.", config.max_tries)

    return Packing(grains=_rebuild_grains(), box=box)
