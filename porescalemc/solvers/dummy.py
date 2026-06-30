"""Convenience dummy solver for examples and tests."""

from __future__ import annotations

from porescalemc.solvers.packing import PackingStatsSolver


class DummySolver(PackingStatsSolver):
    """Alias of :class:`PackingStatsSolver` with a deliberately simple name."""

