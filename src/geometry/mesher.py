"""Gmsh-based mesh generation for periodic porous-media geometries.

Uses the Gmsh Python API with OpenCASCADE (OCC) kernel so that Boolean
operations (fragment, fuse, cut) give a conformal, topologically clean mesh.

The module is a *soft dependency*: if ``gmsh`` is not installed the import
succeeds and all public classes raise ``ImportError`` at construction time.
This keeps the rest of the package import-clean on headless/minimal installs.

Usage example
-------------
    from porescalemc.geometry.mesher import GmshPackingMesher
    from porescalemc.geometry.grains import Packing, Grain

    packing = Packing(box=(1.0, 1.0, 1.0), grains=[Grain.sphere(center=[0.5, 0.5, 0.5], radius=0.2)])
    mesher = GmshPackingMesher(mesh_size=0.05, periodic=True)
    mesher.build(packing)
    mesher.write("packing.msh")
    mesher.finalize()
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Sequence

import numpy as np

from porescalemc.geometry.grains import Grain, Packing

_log = logging.getLogger(__name__)

try:
    import gmsh  # type: ignore
    _GMSH_AVAILABLE = True
except ImportError:
    _GMSH_AVAILABLE = False

# Physical-group names written into the mesh for downstream solvers.
_FLUID_TAG    = "fluid"
_SOLID_TAG    = "solid"
_PERIODIC_TAG = "periodic"


class GmshPackingMesher:
    """Generate a conformal 3-D mesh from a :class:`~porescalemc.geometry.grains.Packing`.

    Parameters
    ----------
    mesh_size : float
        Target element size (same units as the packing box).
    periodic : bool
        If True, enforce face-pairing periodic boundary conditions so that
        the mesh can be used with a periodic Stokes/diffusion solver.
    order : int
        Element order (1 = linear tet, 2 = quadratic tet).
    algo_3d : int
        Gmsh 3-D meshing algorithm tag.  4 = Frontal-Delaunay (robust for
        highly porous geometries); 1 = Delaunay.
    verbosity : int
        Gmsh verbosity level (0 = silent, 5 = very verbose).
    """

    def __init__(
        self,
        mesh_size: float = 0.05,
        periodic: bool = True,
        order: int = 1,
        algo_3d: int = 4,
        verbosity: int = 0,
    ):
        if not _GMSH_AVAILABLE:
            raise ImportError(
                "gmsh is not installed. Install it with: pip install gmsh"
            )
        self.mesh_size = mesh_size
        self.periodic = periodic
        self.order = order
        self.algo_3d = algo_3d
        self.verbosity = verbosity
        self._initialized = False

    # ------------------------------------------------------------------
    # Public interface

    def build(self, packing: Packing) -> None:
        """Build the OCC geometry and mesh for *packing*.

        This initialises Gmsh, creates the box domain and all grain
        volumes, fragments them (conformal intersections), assigns physical
        groups, applies periodicity, and runs the mesher.  Call
        :meth:`write` to save the result and :meth:`finalize` to release
        resources.
        """
        self._init_gmsh()
        Lx, Ly, Lz = packing.box

        # --- Bounding box ---
        box_tag = gmsh.model.occ.addBox(0.0, 0.0, 0.0, Lx, Ly, Lz)

        # --- Grains (spheres / axis-aligned ellipsoids) ---
        grain_tags: list[int] = []
        for grain in packing.grains:
            cx, cy, cz = grain.center
            rx, ry, rz = grain.radii
            if rx == ry == rz:
                tag = gmsh.model.occ.addSphere(cx, cy, cz, rx)
            else:
                # Sphere + dilate gives an axis-aligned ellipsoid.
                tag = gmsh.model.occ.addSphere(cx, cy, cz, 1.0)
                gmsh.model.occ.dilate([(3, tag)], cx, cy, cz, rx, ry, rz)
            grain_tags.append(tag)

        # --- Clip grains to the box (intersection, not subtraction) ---
        # fragment(toolDimTags, objectDimTags) makes a conformal split.
        grain_pairs = [(3, t) for t in grain_tags]
        box_pair    = [(3, box_tag)]

        out_map: list[list[tuple[int, int]]]
        out_vols, _ = gmsh.model.occ.fragment(box_pair, grain_pairs)
        gmsh.model.occ.synchronize()

        # --- Identify fluid and solid volumes ---
        all_vols = [tag for dim, tag in gmsh.model.getEntities(3)]
        fluid_vols, solid_vols = self._classify_volumes(all_vols, packing)

        # --- Physical groups ---
        if fluid_vols:
            gmsh.model.addPhysicalGroup(3, fluid_vols, name=_FLUID_TAG)
        if solid_vols:
            gmsh.model.addPhysicalGroup(3, solid_vols, name=_SOLID_TAG)
        self._add_boundary_physical_groups(Lx, Ly, Lz)

        # --- Periodicity ---
        if self.periodic:
            self._enforce_periodicity(Lx, Ly, Lz)

        # --- Mesh size and generation ---
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", self.mesh_size * 0.5)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", self.mesh_size)
        gmsh.option.setNumber("Mesh.Algorithm3D", self.algo_3d)
        gmsh.model.mesh.generate(3)
        if self.order > 1:
            gmsh.model.mesh.setOrder(self.order)
        _log.info(
            "Meshing complete: %d fluid vols, %d solid vols",
            len(fluid_vols), len(solid_vols),
        )

    def write(self, path: str | Path) -> None:
        """Write mesh to *path* (format inferred from extension, e.g. ``.msh``)."""
        gmsh.write(str(path))
        _log.info("Mesh written to %s", path)

    def finalize(self) -> None:
        """Release Gmsh resources.  Must be called after :meth:`build`/:meth:`write`."""
        if self._initialized:
            gmsh.finalize()
            self._initialized = False

    # ------------------------------------------------------------------
    # Internal helpers

    def _init_gmsh(self) -> None:
        if self._initialized:
            gmsh.finalize()
        gmsh.initialize()
        gmsh.option.setNumber("General.Verbosity", self.verbosity)
        gmsh.model.add("packing")
        self._initialized = True

    def _classify_volumes(
        self, all_vols: list[int], packing: Packing
    ) -> tuple[list[int], list[int]]:
        """Sort volume tags into fluid and solid by centroid containment."""
        fluid_vols: list[int] = []
        solid_vols: list[int] = []
        for vtag in all_vols:
            cx, cy, cz = gmsh.model.occ.getCenterOfMass(3, vtag)
            if _point_in_any_grain(cx, cy, cz, packing.grains):
                solid_vols.append(vtag)
            else:
                fluid_vols.append(vtag)
        return fluid_vols, solid_vols

    def _add_boundary_physical_groups(self, Lx: float, Ly: float, Lz: float) -> None:
        """Add physical surface groups for the six box faces."""
        face_names = {
            "x_minus": (1, 0.0,  0.0,  0.0,  1.0),
            "x_plus":  (1, Lx,   0.0,  0.0,  1.0),
            "y_minus": (2, 0.0,  0.0,  0.0,  1.0),
            "y_plus":  (2, 0.0,  Ly,   0.0,  1.0),
            "z_minus": (3, 0.0,  0.0,  0.0,  1.0),
            "z_plus":  (3, 0.0,  0.0,  Lz,   1.0),
        }
        # Simpler approach: find surfaces by normal direction using BoundingBox
        surfs = [tag for dim, tag in gmsh.model.getEntities(2)]
        tol = min(Lx, Ly, Lz) * 1e-3
        face_groups: dict[str, list[int]] = {
            "x_minus": [], "x_plus": [], "y_minus": [], "y_plus": [],
            "z_minus": [], "z_plus": [],
        }
        for stag in surfs:
            xm, ym, zm, xM, yM, zM = gmsh.model.getBoundingBox(2, stag)
            # A flat box face spans two full dimensions and has no extent in the third
            if abs(xM - xm) < tol and abs(yM - ym) > tol and abs(zM - zm) > tol:
                if xm < tol:
                    face_groups["x_minus"].append(stag)
                elif abs(xm - Lx) < tol:
                    face_groups["x_plus"].append(stag)
            elif abs(yM - ym) < tol and abs(xM - xm) > tol and abs(zM - zm) > tol:
                if ym < tol:
                    face_groups["y_minus"].append(stag)
                elif abs(ym - Ly) < tol:
                    face_groups["y_plus"].append(stag)
            elif abs(zM - zm) < tol and abs(xM - xm) > tol and abs(yM - ym) > tol:
                if zm < tol:
                    face_groups["z_minus"].append(stag)
                elif abs(zm - Lz) < tol:
                    face_groups["z_plus"].append(stag)
        for name, tags in face_groups.items():
            if tags:
                gmsh.model.addPhysicalGroup(2, tags, name=name)

    def _enforce_periodicity(self, Lx: float, Ly: float, Lz: float) -> None:
        """Set mesh periodicity (master-slave face pairing) in all three directions."""
        # Translation vectors for each axis pair.
        translations = [
            ([Lx, 0.0, 0.0], "x_minus", "x_plus"),
            ([0.0, Ly, 0.0], "y_minus", "y_plus"),
            ([0.0, 0.0, Lz], "z_minus", "z_plus"),
        ]
        for tvec, slave_name, master_name in translations:
            slave_tags  = _physical_group_tags(2, slave_name)
            master_tags = _physical_group_tags(2, master_name)
            if slave_tags and master_tags:
                # Affine matrix: 4×4 row-major, last row = [0,0,0,1]
                affine = [
                    1.0, 0.0, 0.0, tvec[0],
                    0.0, 1.0, 0.0, tvec[1],
                    0.0, 0.0, 1.0, tvec[2],
                    0.0, 0.0, 0.0, 1.0,
                ]
                for stag in slave_tags:
                    gmsh.model.mesh.setPeriodic(2, [stag], master_tags, affine)


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------

def _point_in_any_grain(
    cx: float, cy: float, cz: float, grains: Sequence[Grain]
) -> bool:
    """Return True if (cx, cy, cz) is inside any grain ellipsoid."""
    for g in grains:
        gx, gy, gz = g.center
        rx, ry, rz = g.radii
        val = (
            ((cx - gx) / rx) ** 2
            + ((cy - gy) / ry) ** 2
            + ((cz - gz) / rz) ** 2
        )
        if val <= 1.0:
            return True
    return False


def _physical_group_tags(dim: int, name: str) -> list[int]:
    """Return entity tags in a physical group by name, or [] if not found."""
    try:
        tag = gmsh.model.getPhysicalGroupsForName(name)
        if not tag:
            return []
        # tag is a list of (dim, group_tag) tuples matching the name
        for d, gt in tag:
            if d == dim:
                return list(gmsh.model.getEntitiesForPhysicalGroup(d, gt))
    except Exception:
        pass
    return []
