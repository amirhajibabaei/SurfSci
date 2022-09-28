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
    def are_equal(x, y):
        z1, d1 = x
        z2, d2 = y
        n1 = z1.shape[0]
        n2 = z2.shape[0]
        n = min([n1, n2])
        if n1 != n2:
            # This could happen due to numerical errors
            # but the difference shouldn't be more than
            # one particle.
            assert abs(n1 - n2) < 2
        return np.allclose(z1[:n], z2[:n]) and np.allclose(d1[:n], d2[:n])

    if (atoms1.numbers != atoms2.numbers).any():
        return False
    else:
        a, nl = get_nntypes(atoms1, cutoff=cutoff)
        b, _ = get_nntypes(atoms2, nl=nl)
        res = all([are_equal(x, y) for x, y in zip(a, b)])
        return res
