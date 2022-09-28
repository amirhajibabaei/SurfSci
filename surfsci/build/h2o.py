# +
"""
Functions for generating H2O structures.
All the functions that return an Atoms object
should sort atoms such that an Oxygen (O) is
followed by the two Hydrogen (H) atoms which
are bonded with the Oxygen:
   O
   H   molecule 1
   H

   O
   H   molecule 2
   H
   ...

"""
import io
from ase.io import read
from ase import Atoms
from ase.neighborlist import neighbor_list
import numpy as np


def tetragonal_ice_unitcell() -> Atoms:
    return read(_tetragonal_ice_xyz, format="extxyz")


def block_of_ice(
    lxlylz: tuple[float, float, float],
    ciel: tuple[bool, bool, bool] = (False, False, False),
):
    """
    It will generate a block of ice using the tetragonal
    unit cell.

    Args:
        lxlylz:  target block lengths
        ciel:    if True, it is allowed to surpass the
                 length by a maximum of one unit cell
    """
    uc = tetragonal_ice_unitcell()
    size = [
        int(np.ceil(l / u) if c else np.floor(l / u))
        for u, l, c in zip(uc.cell.diagonal(), np.asarray(lxlylz), ciel)
    ]
    return uc.repeat(size)


def unwrap_oh_bounds(atoms: Atoms, bmax: float = 1.2) -> None:
    # TODO: better algorithm
    nl = neighbor_list("ijd", atoms, cutoff=bmax)
    for i, j, d in zip(*nl):
        a = atoms[i]
        b = atoms[j]
        if a.symbol == "O" and b.symbol == "H":
            r_mic = atoms.get_distance(i, j, mic=True, vector=True)
            b.position = a.position + r_mic


def append_h2o(atoms: Atoms, length: float, vacuum: float, ciel: bool = True) -> Atoms:
    """
    Append a block of H2O molecules (tetragonal ice)
    along the z axis of the given structure.

    Args:
        atoms:    input atoms (often a surface)
        length:   target length of the ice block along z
        vacuum:   vacuum on both sides of the ice block
        ciel:     if True, size is allowed to exceed
                  "length" by one unitcell of ice
    """
    par = atoms.cell.cellpar()
    assert np.allclose(par[3:], 90)

    # create H2O
    lx, ly = par[:2]
    h2o = block_of_ice((lx, ly, length), ciel=(False, False, ciel))
    h2o.cell[0, 0] = lx
    h2o.cell[1, 1] = ly
    h2o.center(axis=(0, 1))

    #
    pos1 = atoms.positions
    pos1 -= pos1.min(axis=0)
    pos2 = h2o.positions
    pos2[:, 2] += -pos2[:, 2].min() + pos1[:, 2].max() + vacuum
    lz = pos2[:, 2].max() + vacuum

    result = Atoms(
        numbers=np.concatenate([atoms.numbers, h2o.numbers]),
        positions=np.concatenate([pos1, pos2]),
        cell=[lx, ly, lz],
        pbc=atoms.pbc,
    )
    return result


_tetragonal_ice_xyz = io.StringIO(
    """6
Lattice="3.413821 0.0 0.0 0.0 3.413821 0.0 0.0 0.0 3.381228" Properties=species:S:1:pos:R:3 pbc="T T T"
O        2.56036575       2.56036575       2.83648174      -0.00000000
H        2.00371834       2.00371834       2.23536026       0.00000000
H        3.11701316       3.11701316       2.23536026       0.00000000
O        0.85345525       0.85345525       1.14586774      -0.00000000
H        0.29680784       1.41010266       0.54474626       0.00000000
H        1.41010266       0.29680784       0.54474626       0.00000000
"""  # noqa: E501
)
