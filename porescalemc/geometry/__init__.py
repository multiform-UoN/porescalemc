"""Geometry sub-package: grain data structures, placement, and field generation."""

from __future__ import annotations

from porescalemc.geometry.grains import Grain, Packing
from porescalemc.geometry.placement import sample_packing, sample_packing_from_name

# GmshPackingMesher is a soft dependency — available only when gmsh is installed.
try:
    from porescalemc.geometry.mesher import GmshPackingMesher
    __all__ = ["Grain", "Packing", "sample_packing", "sample_packing_from_name", "GmshPackingMesher"]
except ImportError:
    __all__ = ["Grain", "Packing", "sample_packing", "sample_packing_from_name"]
