# +
from __future__ import annotations

from ase import Atoms
from ase.build import surface
from ase.spacegroup import crystal

from .transform import right_angle_cell
from .h2o import append_h2o


def nacl_unitcell(a: float | None = None) -> Atoms:
    """
    Creates a NaCl unit cell with spacegroup 225.

    Args:
        a:   optional lattice spacing (if None, a = 5.6402)

    """
    if a is None:
        a = 5.6402
    unitcell = crystal(
        symbols=("Na", "Cl"),
        basis=((0, 0, 0), (0.5, 0.5, 0.5)),
        spacegroup=225,
        cellpar=(a, a, a, 90, 90, 90),
    )
    return unitcell


def nacl_surface(
    a: float | None = None,
    indices: tuple[int, int, int] = (1, 1, 1),
    size: tuple[int, int, int] = (1, 1, 1),
    right_angle: bool = True,
    vacuum: float | None = None,
    h2o: float | None = None,
    periodic: bool = True,
) -> Atoms:
    """
    Creates NaCl surface (based on spacegroup 225)
    where z-axis will be normal to the surface.

    Args:
        a:            lattice spacing
        indices:      miller indices
        size:         (nx, ny, nz) where nz is the number
                      of layers
        right_angle:  if True, ensures a cell with right angles
        vacuum:       optional vacuum on top of surface or between
                      the surface and water (if finite h2o)
        h2o:          if given a block of ice will be deposited
                      on the surface
        periodic:     periodic boundary condition (along z)

    Restrictions:
        If finite h2o is dpecified, finite vacuum should
        also be specified.

    """
    atoms = surface(
        lattice=nacl_unitcell(a),
        indices=indices,
        layers=size[2],
        periodic=periodic,
        vacuum=vacuum,
    ).repeat((*size[:2], 1))
    if right_angle:
        atoms = right_angle_cell(atoms)
    if h2o is not None:
        assert vacuum is not None
        atoms = append_h2o(atoms, h2o, vacuum)
    return atoms
