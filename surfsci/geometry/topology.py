# +
from itertools import combinations
from typing import Sequence

import numpy as np
from ase import Atoms
from ase.neighborlist import NeighborList

from ..assertions import in_range_rtol
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
                    if atoms[j].symbol == bb and in_range_rtol(d, bd, rtol):
                        bonds_a.append((t, j))
        bonds.append(bonds_a)
    return bonds


def find_angles(
    atoms: Atoms,
    bonds: list[list[tuple[int, int]]],
    angle_types: Sequence[tuple[int, int, float]],
    rtol: float = 1.1,
):
    """
    Args:
        atoms:        ase.Atoms object
        bonds:        output of find_bonds
        angle_types:  a sequence of angle types e.g.
                        [(0, 0, 1.), ...]
                      where the first two are indices of two bond
                      types (which is used for creating the bonds)
                      and the third is a float indicating the angle.
        rtol:         relative tolerance for angle

    Returns:
        a list of angles for each atom.
        [
            [(t1, j1, k1), ...],   <- angle j1--0--k1
            [(t2, j2, k2), ...],   <- bonds j2--1--k2
                 ...                      etc.
        ]
        where an angle is encoded by tuple of (t, j, k) in which
        t is the type of angle (index in angle_types) and
        j,k are the indices of two other atoms forming the angle.

    """
    angles = []
    for a, bonds_a in enumerate(bonds):
        angles_a = []
        for b1, b2 in combinations(bonds_a, 2):
            type_ = {b1[0], b2[0]}  # set
            for t, (bt1, bt2, angle) in enumerate(angle_types):
                if {bt1, bt2} == type_:  # set equality
                    angle_ = atoms.get_angle(b1[1], a, b2[1])
                    if in_range_rtol(angle, angle_, rtol):
                        angles_a.append((t, b1[1], b2[1]))
        angles.append(angles_a)
    return angles
