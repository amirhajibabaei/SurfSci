# +
import numpy as np
from ase import Atoms
from ase.neighborlist import NeighborList


def get_neighbors(
    atoms: Atoms, nl: NeighborList, i: int
) -> tuple[np.ndarray, np.ndarray]:
    j, off = nl.get_neighbors(i)
    r = (
        atoms.positions[j]
        - atoms.positions[i]
        + (off[..., None] * atoms.cell).sum(axis=1)
    )
    return j, r
