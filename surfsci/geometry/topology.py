# +
from typing import Sequence

import numpy as np
from ase import Atoms
from ase.neighborlist import NeighborList

from .nlutil import get_neighbors


def find_bonds(
    atoms: Atoms, bond_types: Sequence[tuple[str, str, float]], rtol: float = 0.1
) -> list[list[tuple[int, int]]]:
    """
    Args:
        atoms:        ase.Atoms object
        bond_types:   a sequence of bond types e.g.
                        [("O", "H", 1.), ...]
                      where the first two are atomic symbols
                      and the third is a float indicating distance.
        rtol:         relative tolerance for distance

    Returns:
        a list of bonds for each atom.
        [
            [(t1, j1), ...],   <- bonds connected to atom 0
            [(t2, j2), ...],   <- bonds connected to atom 1
                 ...                      etc.
        ]
        where a bond is encoded by tuple of (t, j) in which
        t is the type of bond (index in bond_types) and
        j is the index of the other atom in the bond.

    """
    cutoff = (1 + rtol) * max([d for _, _, d in bond_types])
    nl = NeighborList(
        len(atoms) * [cutoff / 2],
        skin=0.0,
        # sorted=True,
        self_interaction=False,
        bothways=True,
    )
    nl.update(atoms)
    bonds = []
    for a in atoms:
        ind, r = get_neighbors(atoms, nl, a.index)
        dis = np.linalg.norm(r, axis=1)
        nei = zip(ind, dis)
        bonds_a = []
        for t, (ba, bb, bd) in enumerate(bond_types):
            if a.symbol == ba:
                for j, d in nei:
                    if (
                        atoms[j].symbol == bb
                        and d > bd * (1 - rtol)
                        and d < bd * (1 + rtol)
                    ):
                        bonds_a.append((t, j))
        bonds.append(bonds_a)
    return bonds
