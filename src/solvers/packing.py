"""Algebraic packing statistics as a pure-Python solver.

These quantities are useful as dummy QoIs for MLMC development and are also
standard low-order descriptors of random media: porosity, number density,
specific surface area, size-distribution moments, and nearest-neighbour
statistics.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from porescalemc.geometry.grains import Grain, Packing, ellipsoid_distance


PACKING_QOI_NAMES = (
    "porosity",
    "solid_fraction",
    "number_density",
    "specific_surface_area",
    "mean_equivalent_radius",
    "radius_cv",
    "mean_nearest_center_distance",
    "mean_nearest_normalized_gap",
)


def _equivalent_radius(grain: Grain) -> float:
    """Radius of the sphere with the same volume as ``grain``."""
    return float((3.0 * grain.volume / (4.0 * math.pi)) ** (1.0 / 3.0))


def _surface_area(grain: Grain) -> float:
    """Grain surface area.

    Spheres are exact.  Ellipsoids use the Knud Thomsen approximation::

        S ≈ 4π · ((ab)^p + (ac)^p + (bc)^p) / 3)^(1/p),   p = 1.6075

    This approximation has error < 1.06% for axis ratios up to 10:1,
    sufficient for diagnostic/statistical QoIs; no special functions needed.
    """
    if not grain.is_ellipsoid:
        r = float(grain.transform[0])
        return 4.0 * math.pi * r**2

    a, b, c = grain.radii
    p = 1.6075  # Knud Thomsen exponent
    return float(
        4.0
        * math.pi
        * (((a * b) ** p + (a * c) ** p + (b * c) ** p) / 3.0) ** (1.0 / p)
    )


def _nearest_statistics(packing: Packing) -> tuple[float, float]:
    """Return mean nearest centre distance and mean nearest normalised gap.

    Note: this is O(N²) in the number of grains.  For diagnostic QoIs on
    typical packings (N < 200) the cost is negligible; for N > 1000 replace
    with scipy.spatial.cKDTree for O(N log N) performance.
    """
    n = packing.n_grains
    if n < 2:
        return math.nan, math.nan

    nearest_center = []
    nearest_gap = []
    box = packing.box
    for i, gi in enumerate(packing.grains):
        center_dists = []
        normalized_gaps = []
        for j, gj in enumerate(packing.grains):
            if i == j:
                continue
            dvec = gj.center - gi.center
            dist = float(np.linalg.norm(dvec))
            center_dists.append(dist)
            # Normalised gap = ellipsoid_distance - 1; negative means overlap.
            normalized_gaps.append(ellipsoid_distance(gi, gj, box=box) - 1.0)
        nearest_center.append(min(center_dists))
        nearest_gap.append(min(normalized_gaps))

    return float(np.mean(nearest_center)), float(np.mean(nearest_gap))


def packing_statistics(packing: Packing) -> dict[str, float]:
    """Compute named algebraic statistics for ``packing``."""
    radii = np.array([_equivalent_radius(g) for g in packing.grains], dtype=float)
    if radii.size == 0:
        mean_radius = math.nan
        radius_cv = math.nan
    else:
        mean_radius = float(np.mean(radii))
        radius_cv = float(np.std(radii, ddof=0) / mean_radius) if mean_radius > 0 else math.nan

    surface_area = sum(_surface_area(g) for g in packing.grains)
    mean_nn_dist, mean_nn_gap = _nearest_statistics(packing)
    solid_fraction = 1.0 - packing.porosity

    return {
        "porosity": float(packing.porosity),
        "solid_fraction": float(solid_fraction),
        "number_density": float(packing.n_grains / packing.box_volume),
        "specific_surface_area": float(surface_area / packing.box_volume),
        "mean_equivalent_radius": mean_radius,
        "radius_cv": radius_cv,
        "mean_nearest_center_distance": mean_nn_dist,
        "mean_nearest_normalized_gap": mean_nn_gap,
    }


@dataclass
class PackingStatsSolver:
    """Solver adapter that returns algebraic packing statistics."""

    packing: Packing | None = None
    qoi_names: tuple[str, ...] = PACKING_QOI_NAMES

    def setup(self, packing: Packing) -> None:
        self.packing = packing

    def solve(self) -> np.ndarray:
        if self.packing is None:
            raise RuntimeError("PackingStatsSolver.setup must be called before solve.")
        stats = packing_statistics(self.packing)
        return np.array([stats[name] for name in self.qoi_names], dtype=float)

    def close(self) -> None:
        return None
