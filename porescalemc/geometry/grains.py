"""
Grain and Packing data structures, plus geometry primitives.

A ``Grain`` is either a sphere or an oriented ellipsoid:

- Sphere:    ``transform`` is a 1-element array ``[r]`` (radius).
- Ellipsoid: ``transform`` is a 9-element array storing the flattened 3×3
             transformation matrix M that maps the unit ball to the ellipsoid::

                 ellipsoid = { x0 + M·u  :  |u| ≤ 1 }

  The semi-axis lengths are the column norms of M; orientation is encoded in
  M's columns.  This is the same convention as the original geomodule.py.

Distance metric
---------------
The normalised ellipsoid distance between two grains g1, g2 is::

    d(g1, g2) = |c2 - c1| / (s(g1, n̂) + s(g2, n̂))

where n̂ = (c2 - c1)/|c2 - c1| and s(g, n̂) is the stretching length of
grain g in direction n̂ (i.e., the radius of g in that direction).

d < 1  → grains overlap
d = 1  → grains just touch
d > 1  → grains are separated
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np


# ---------------------------------------------------------------------------
# Grain
# ---------------------------------------------------------------------------

@dataclass
class Grain:
    """A single grain: sphere or oriented ellipsoid.

    Parameters
    ----------
    center : np.ndarray, shape (3,)
        Position of the grain centre.
    transform : np.ndarray
        Shape (1,) for a sphere (scalar radius ``r``), or shape (9,) for an
        ellipsoid (flattened 3×3 transformation matrix M, column-major order).
    """

    center: np.ndarray    # shape (3,)
    transform: np.ndarray  # shape (1,) or (9,)

    # ------------------------------------------------------------------ #
    # Construction helpers                                                 #
    # ------------------------------------------------------------------ #

    @classmethod
    def sphere(cls, center: np.ndarray | list, radius: float) -> Grain:
        """Create a sphere grain."""
        return cls(
            center=np.asarray(center, dtype=float),
            transform=np.array([float(radius)]),
        )

    @classmethod
    def ellipsoid(cls, center: np.ndarray | list, M: np.ndarray) -> Grain:
        """Create an ellipsoid grain from a 3×3 transformation matrix M."""
        return cls(
            center=np.asarray(center, dtype=float),
            transform=np.asarray(M, dtype=float).ravel(),
        )

    # ------------------------------------------------------------------ #
    # Properties                                                           #
    # ------------------------------------------------------------------ #

    @property
    def is_ellipsoid(self) -> bool:
        """True if this grain is an oriented ellipsoid, False if sphere."""
        return len(self.transform) == 9

    @property
    def radii(self) -> np.ndarray:
        """Semi-axis lengths as a (3,) array.

        For a sphere, all three values equal the radius.
        For an ellipsoid, they are the column norms of M.
        """
        if self.is_ellipsoid:
            M = self.transform.reshape(3, 3)
            # Column norms of M give the semi-axis lengths (M maps unit sphere to ellipsoid).
            return np.sqrt((M ** 2).sum(axis=0))
        return np.full(3, self.transform[0])

    @property
    def max_radius(self) -> float:
        """Maximum semi-axis length."""
        return float(self.radii.max())

    @property
    def volume(self) -> float:
        """Exact grain volume."""
        if self.is_ellipsoid:
            M = self.transform.reshape(3, 3)
            # Volume = (4π/3)|det M|; det M = product of semi-axes for aligned ellipsoids.
            return 4.0 / 3.0 * math.pi * abs(float(np.linalg.det(M)))
        return 4.0 / 3.0 * math.pi * float(self.transform[0]) ** 3

    def as_array(self) -> np.ndarray:
        """Concatenated ``[center, transform]`` — legacy-compatible representation."""
        return np.concatenate([self.center, self.transform])

    def transformation_matrix(self) -> np.ndarray:
        """3×3 matrix M (identity*r for spheres, full M for ellipsoids)."""
        if self.is_ellipsoid:
            return self.transform.reshape(3, 3)
        r = float(self.transform[0])
        return np.eye(3) * r


# ---------------------------------------------------------------------------
# Packing
# ---------------------------------------------------------------------------

@dataclass
class Packing:
    """A collection of grains inside a rectangular periodic box.

    Parameters
    ----------
    grains : list[Grain]
        All grains in the domain.
    box : np.ndarray, shape (3,)
        Domain lengths [Lx, Ly, Lz].  The box is centred at the origin, so
        grain centres live in [-Lx/2, Lx/2] × [-Ly/2, Ly/2] × [-Lz/2, Lz/2].
    """

    grains: list[Grain]
    box: np.ndarray  # shape (3,)

    @property
    def n_grains(self) -> int:
        """Number of grains."""
        return len(self.grains)

    @property
    def box_volume(self) -> float:
        """Total domain volume."""
        return float(np.prod(self.box))

    @property
    def solid_volume(self) -> float:
        """Sum of grain volumes (may exceed box_volume for overlapping grains)."""
        return sum(g.volume for g in self.grains)

    @property
    def porosity(self) -> float:
        """Estimated porosity = 1 - solid_volume / box_volume.

        This is an estimate: it ignores grain–grain and grain–boundary overlaps.
        Use the Fourier field for an accurate value.
        """
        return 1.0 - min(1.0, self.solid_volume / self.box_volume)

    @property
    def centers(self) -> np.ndarray:
        """All grain centres as an (N, 3) array."""
        if not self.grains:
            return np.empty((0, 3))
        return np.array([g.center for g in self.grains])

    def __len__(self) -> int:
        return self.n_grains


# ---------------------------------------------------------------------------
# Geometry primitives
# ---------------------------------------------------------------------------

def stretching_length(grain: Grain, direction: np.ndarray) -> float:
    """Radius of ``grain`` along ``direction`` (must be a unit vector).

    For a sphere: returns the scalar radius.
    For an ellipsoid with matrix M: returns ||M^T · n̂||, which is the
    maximum extent of the ellipsoid in direction n̂.

    Parameters
    ----------
    grain : Grain
    direction : np.ndarray, shape (3,)
        Unit direction vector.
    """
    n = np.asarray(direction, dtype=float)
    if grain.is_ellipsoid:
        M = grain.transform.reshape(3, 3)
        # ||M^T n|| is the support function of the ellipsoid: it gives the half-width
        # of the ellipsoid in direction n (i.e., its "radius" in that direction).
        return float(np.linalg.norm(M.T @ n))
    return float(grain.transform[0])


def _periodic_displacement(c1: np.ndarray, c2: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Minimum-image displacement vector c2 - c1 under periodic BCs."""
    d = c2 - c1
    # np.round maps each component to the nearest integer; dividing by box beforehand
    # then multiplying back gives the closest periodic image of c2 relative to c1.
    return d - box * np.round(d / box)


def ellipsoid_distance(
    g1: Grain,
    g2: Grain,
    box: np.ndarray | None = None,
) -> float:
    """Normalised centre-to-surface distance between two grains.

    Returns
    -------
    float
        d < 1 → overlapping, d = 1 → touching, d > 1 → separated.
        Returns np.inf if grains are coincident (avoid division by zero).

    Parameters
    ----------
    g1, g2 : Grain
    box : np.ndarray, shape (3,), optional
        If given, apply minimum-image convention for periodic BCs.
    """
    if box is not None:
        d = _periodic_displacement(g1.center, g2.center, np.asarray(box))
    else:
        d = g2.center - g1.center

    dist = float(np.linalg.norm(d))
    if dist < 1e-15:
        return np.inf

    # Unit vector pointing from g1 to g2 (direction of closest approach).
    n = d / dist
    # Sum of stretching lengths = total "diameter" of the two grains along n.
    # Normalising the distance by this sum gives the dimensionless overlap metric.
    denom = stretching_length(g1, n) + stretching_length(g2, n)
    if denom < 1e-15:
        return np.inf
    return dist / denom


def sample_random_orientation(radii: list[float], rng: np.random.Generator) -> np.ndarray:
    """Sample a random orientation matrix M using the Haar measure.

    Generates a uniformly distributed random rotation Q (from the
    Gram-Schmidt + sign-normalisation method of Maris Ozols 2009), then
    scales each column by the corresponding radius::

        M = Q · diag(radii)

    This produces a 3×3 matrix whose column norms equal ``radii`` and whose
    orientation is uniformly distributed on SO(3).

    Parameters
    ----------
    radii : list of 3 floats
        Semi-axis lengths.
    rng : numpy.random.Generator

    Returns
    -------
    M_flat : np.ndarray, shape (9,)
        Flattened transformation matrix.
    """
    # Fill a 3×3 matrix with i.i.d. standard normals; QR decomposition then
    # gives a Haar-distributed orthogonal matrix Q.
    A = rng.standard_normal((3, 3))
    Q, R = np.linalg.qr(A)
    # The sign of each column of Q is ambiguous; multiplying by sign(diag(R))
    # makes the representation unique and preserves Haar measure (Ozols 2009).
    Q = Q @ np.diag(np.sign(np.diag(R)))
    # Force det = +1 so Q ∈ SO(3); for an ellipsoid the sign doesn't matter
    # geometrically (centrosymmetric), but it keeps the math consistent.
    Q = Q / np.linalg.det(Q)
    M = Q @ np.diag(np.asarray(radii, dtype=float))
    return M.ravel()


def grains_to_array(grains: list[Grain]) -> np.ndarray:
    """Stack grains as (N, D) array where D = 4 (sphere) or 12 (ellipsoid)."""
    if not grains:
        return np.empty((0, 4))
    return np.array([g.as_array() for g in grains])
