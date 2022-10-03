# +
from __future__ import annotations

from ase import Atoms
from ase.build import surface as ase_surface, add_vacuum

from .transform import right_angle_cell, conventional_cell
from .h2o import append_h2o


def surface(
    unitcell: Atoms,
    indices: tuple[int, int, int] = (0, 0, 1),
    size: tuple[int, int, int] = (1, 1, 1),
    conventional: bool = True,
    right_angle: bool = True,
    vacuum: float | None = None,
    h2o: float | None = None,
    periodic: bool = True,
) -> Atoms:
    """
    Creates a surface along the given miller indices
    where z-axis will be normal to the surface.
    Optionally it can add H2O on top of the surface.

    Args:
        unitcell:     the unit-cell Atoms object
        indices:      miller indices of the surface normal
                      to the z axis
        size:         (nx, ny, nz) where nz is the number
                      of layers
        conventional: if True, uses the conventional unit cell
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
    if conventional:
        unitcell = conventional_cell(unitcell)

    atoms = ase_surface(
        lattice=unitcell,
        indices=indices,
        layers=size[2],
        periodic=periodic,
        vacuum=None,
    ).repeat((*size[:2], 1))

    if right_angle:
        atoms = right_angle_cell(atoms)

    if h2o is not None:
        assert vacuum is not None
        atoms = append_h2o(atoms, h2o, vacuum)
    elif vacuum is not None:
        add_vacuum(atoms, vacuum)

    return atoms
