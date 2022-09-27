# +
"""
Functions for testing equality of structures.

"""
import numpy as np
from ase import Atoms
from ase.neighborlist import NeighborList


def get_nntypes(
    atoms: Atoms, cutoff: float = 5.0, nl: NeighborList = None
) -> tuple[list[tuple[np.ndarray, np.ndarray]], NeighborList]:
    """
    For each atom, finds atomic types of its
    neighbors sorted based on their distance.
    A list of (types, distances) is returned
    along with generated/updated neighborlist.

    """
    n = len(atoms)
    if nl is None:
        nl = NeighborList(
            n * [cutoff / 2],
            skin=0.0,
            # sorted=True,
            self_interaction=False,
            bothways=True,
        )
    nl.update(atoms)
    types = []
    for i in range(n):
        j, off = nl.get_neighbors(i)
        r = (
            atoms.positions[j]
            - atoms.positions[i]
            + (off[..., None] * atoms.cell).sum(axis=1)
        )
        d = np.linalg.norm(r, axis=1)
        arg = np.argsort(d)
        types.append((atoms.numbers[j[arg]], d[arg]))
    return types, nl


def locally_identical(atoms1: Atoms, atoms2: Atoms, cutoff=5.0):
    if (atoms1.numbers != atoms2.numbers).any():
        return False
    else:
        a, nl = get_nntypes(atoms1, cutoff=cutoff)
        b, _ = get_nntypes(atoms2, nl=nl)
        res = all(
            [(x[0] == y[0]).all() and np.allclose(x[1], y[1]) for x, y in zip(a, b)]
        )
        return res
