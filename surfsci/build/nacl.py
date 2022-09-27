# +
from __future__ import annotations

from ase import Atoms
from ase.build import surface
from ase.spacegroup import crystal


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
    vacuum: float | None = None,
    periodic: bool = True,
) -> Atoms:
    """
    Creates NaCl surface (based on spacegroup 225)
    where z-axis will be normal to the surface.

    Args:
        a:          lattice spacing
        indices:    miller indices
        size:       (nx, ny, nz) where nz is the number
                    of layers
        vacuum:     optional vacuum on top of surface
        periodic:   periodic boundary condition (along z)

    """

    surf = surface(
        lattice=nacl_unitcell(a),
        indices=indices,
        layers=size[2],
        periodic=periodic,
        vacuum=vacuum,
    )
    return surf.repeat((*size[:2], 1))
