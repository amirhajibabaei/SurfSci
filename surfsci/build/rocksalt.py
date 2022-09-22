# +
from functools import partial

from ase.spacegroup import crystal
from ase import Atoms

__all__ = ["rocksalt"]


def rocksalt(
    a: float, symbols: tuple[str, str], size: tuple[int, int, int] = (1, 1, 1)
) -> Atoms:
    atoms = crystal(
        symbols=symbols,
        basis=((0, 0, 0), (0.5, 0.5, 0.5)),
        spacegroup=225,
        cellpar=(a, a, a, 90, 90, 90),
    )
    return atoms.repeat(size)


NaCl = partial(rocksalt, a=5.6402, symbols=("Na", "Cl"))
