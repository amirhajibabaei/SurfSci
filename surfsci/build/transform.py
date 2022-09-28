# +
from math import pi, radians

import numpy as np

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.transformations.transformation_abc import AbstractTransformation
from pymatgen.transformations.standard_transformations import (
    ConventionalCellTransformation,
)

from ase import Atoms

from ..geometry.rotations import AxisRotation, PrismRotation
from ..geometry.equallity import locally_identical


def right_angle_cell(atoms: Atoms) -> Atoms:
    """
    Tries to convert a prism cell into cubic by cutting
    and moving the corners. Note that this will NOT
    preserve the indexwise neighborhood relations.
    That means that (i, j) distance may change
    in the generated structure. But this could still
    be usefull if the local chemistry remains the same
    for all atoms. As certain translational symmetry
    is required, this may not allways be possible.
    Thus an error will be raised if the transformation
    is not possible.

    It returns an Atoms object without modifying the input.

    """
    # i. test if the input is qualified
    assert all(atoms.pbc)
    angles = atoms.cell.cellpar()[3:]
    ang = pi / 2 - radians(angles[2])
    if not np.allclose(angles[:2], 90) and ang > 0:
        raise RuntimeError("Only cells with 90, 90, < 90 angles are supported for now!")
        # TODO: generalize

    # ii. create
    rot1 = PrismRotation(atoms.cell)
    rot2 = AxisRotation((0, 0, 1), ang)
    right = Atoms(
        numbers=atoms.numbers,
        positions=rot2(rot1(atoms.positions)),
        cell=rot2(rot1(atoms.cell)),
        pbc=atoms.pbc,
    )
    assert np.allclose(right.cell[2, :2], 0)
    right.cell[0, 1] = 0
    right.wrap()

    # iii. test
    if not locally_identical(atoms, right):
        raise RuntimeError("Transformation is not possible!")
    return right


def apply_transformation(tr: AbstractTransformation, atoms: Atoms) -> Atoms:
    adaptor = AseAtomsAdaptor()
    struc = adaptor.get_structure(atoms)
    c = tr.apply_transformation(struc)
    return adaptor.get_atoms(c)


def conventional_cell(atoms: Atoms) -> Atoms:
    return apply_transformation(ConventionalCellTransformation(), atoms)
