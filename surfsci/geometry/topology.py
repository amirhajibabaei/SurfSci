# +
from __future__ import annotations

from itertools import combinations
from typing import Sequence

import numpy as np
from ase import Atoms
from ase.neighborlist import NeighborList

from ..assertions import in_range_rtol
from .nlutil import get_neighbors


def find_bonds(
    atoms: Atoms, bond_types: Sequence[tuple[str, str, float]], rtol: float = 0.1
) -> list[list[tuple[int, int, int]]]:
    """
    Args:
        atoms:        ase.Atoms object
        bond_types:   a sequence of bond types e.g.
                        [("O", "H", 1.), ...]
                      where the first two are atomic symbols
                      and the third is a float indicating distance.
        rtol:         relative tolerance for distance

    Returns:
        a list of bonds for each atom:
        [
                ...
            [(t1, i, j1), (t2, i, j2), ...], # bonds on atom i
                ...
        ]
        where a bond is encoded by tuple of (t, i, j) in which
        t is the type of bond (index in bond_types) and i, j
        index of atoms in the bond. The index i is the same
        for all bonds of atom i.

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
                        bonds_a.append((t, a.index, j))
        bonds.append(bonds_a)
    return bonds


def find_angles(
    atoms: Atoms,
    bonds: list[list[tuple[int, int, int]]],
    angle_types: Sequence[tuple[int, int, float]],
    rtol: float = 1.1,
) -> list[list[tuple[int, int, int, int]]]:
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
                ...
            [(t1, j1, i, k1), (t2, j2, i, k2), ...], # angles on atom i
                ...
        ]
        where an angle is encoded by tuple of (t, j, i, k) in which
        t is the type of angle (index in angle_types) and j, k are
        the indices of two other atoms forming j--i--k angle.

    """
    angles = []
    for a, bonds_a in enumerate(bonds):
        angles_a = []
        for b1, b2 in combinations(bonds_a, 2):
            type_ = {b1[0], b2[0]}  # set
            for t, (bt1, bt2, angle) in enumerate(angle_types):
                if {bt1, bt2} == type_:  # set equality
                    angle_ = atoms.get_angle(b1[2], a, b2[2])
                    if in_range_rtol(angle, angle_, rtol):
                        angles_a.append((t, b1[2], a, b2[2]))
        angles.append(angles_a)
    return angles


def find_topology(atoms: Atoms, topology: dict, rtol: float = 0.1) -> dict:
    res: dict[str, list] = {"bond": [], "angle": []}
    if "bond" in topology:
        res["bond"] = find_bonds(atoms, topology["bond"], rtol)
    if "angle" in topology:
        res["angle"] = find_angles(atoms, res["bond"], topology["angle"], rtol)
    return res
