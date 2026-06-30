"""
Porosity-to-permeability closure laws.

All laws take a porosity array ε ∈ (0, 1) and a characteristic grain size d,
and return the dimensionless permeability K (or inverse permeability 1/K).

Convention: ε is the *porosity* (void fraction), so ε = 0 → fully solid,
ε = 1 → fully empty.  This is the opposite of the original code, which used
the solid fraction (1 - porosity) in many places.

Kozeny-Carman
^^^^^^^^^^^^^
K = ε³ / (180 · (1 - ε)² / d²)

Ergun
^^^^^
Extends Kozeny-Carman with an inertial correction::

    1/K = 150(1-ε)²/(ε³ d²) + 1.75 Re (1-ε)/(ε³ d)

Wen-Yu
^^^^^^
Empirical correlation for packed beds::

    K = ε^3.65 / (18 · (1-ε) / d²)
"""

from __future__ import annotations

import numpy as np


def kozeny_carman(
    porosity: np.ndarray | float,
    grain_size: float,
    limiter: float = 1e-3,
) -> np.ndarray:
    """Kozeny-Carman permeability.

    K = ε³ d² / (180 (1 − ε)²)

    Parameters
    ----------
    porosity : array-like
        Void fraction ε ∈ (0, 1).
    grain_size : float
        Characteristic diameter d.
    limiter : float
        Minimum porosity value to avoid division by zero.

    Returns
    -------
    np.ndarray
        Permeability (same units as grain_size²).
    """
    # Clamp ε away from 0 and 1: at ε→0 the numerator → 0 faster than denominator,
    # but at ε→1 the denominator → 0 causing K → ∞ (pure fluid — not physical in a packing).
    eps = np.clip(np.asarray(porosity, dtype=float), limiter, 1.0 - limiter)
    return eps ** 3 * grain_size ** 2 / (180.0 * (1.0 - eps) ** 2)


def ergun(
    porosity: np.ndarray | float,
    grain_size: float,
    reynolds: float = 0.0,
    limiter: float = 1e-3,
) -> np.ndarray:
    """Ergun permeability (Kozeny-Carman + inertial term).

    1/K = 150(1-ε)²/(ε³ d²) + 1.75 Re (1-ε)/(ε³ d)

    Parameters
    ----------
    porosity : array-like
    grain_size : float
    reynolds : float
        Reynolds number Re (0 for Stokes flow → reduces to Kozeny-Carman).
    limiter : float

    Returns
    -------
    np.ndarray
    """
    # limiter prevents K→∞ at ε→1 and 1/K→∞ at ε→0.
    eps = np.clip(np.asarray(porosity, dtype=float), limiter, 1.0 - limiter)
    inv_K = (
        150.0 * (1.0 - eps) ** 2 / (eps ** 3 * grain_size ** 2)  # Kozeny-Carman viscous term
        + 1.75 * reynolds * (1.0 - eps) / (eps ** 3 * grain_size)  # Forchheimer inertial correction
    )
    return 1.0 / np.where(inv_K > 0, inv_K, np.inf)


def wen_yu(
    porosity: np.ndarray | float,
    grain_size: float,
    limiter: float = 1e-3,
) -> np.ndarray:
    """Wen-Yu permeability (empirical packed-bed correlation).

    K = ε^3.65 d² / (18 (1 − ε))

    Parameters
    ----------
    porosity : array-like
    grain_size : float
    limiter : float

    Returns
    -------
    np.ndarray
    """
    # limiter: same reason as Kozeny-Carman — prevent singularities at the domain extremes.
    eps = np.clip(np.asarray(porosity, dtype=float), limiter, 1.0 - limiter)
    return eps ** 3.65 * grain_size ** 2 / (18.0 * (1.0 - eps))


def permeability_field(
    porosity: np.ndarray | float,
    law: str,
    grain_size: float,
    **kwargs,
) -> np.ndarray:
    """Apply a named porosity-permeability law to a field array.

    Parameters
    ----------
    porosity : array-like
        Porosity field (void fraction ε).
    law : str
        'kozeny_carman', 'ergun', or 'wen_yu'.
    grain_size : float
        Characteristic grain diameter.
    **kwargs
        Additional keyword arguments forwarded to the law function
        (e.g., ``reynolds=1.0`` for Ergun).

    Returns
    -------
    np.ndarray
        Permeability field.
    """
    porosity = np.asarray(porosity, dtype=float)
    if law == "kozeny_carman":
        return kozeny_carman(porosity, grain_size, **kwargs)
    elif law == "ergun":
        return ergun(porosity, grain_size, **kwargs)
    elif law == "wen_yu":
        return wen_yu(porosity, grain_size, **kwargs)
    else:
        raise ValueError(
            f"Unknown law: {law!r}. Choose 'kozeny_carman', 'ergun', or 'wen_yu'."
        )
