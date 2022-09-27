# +
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.transformations.transformation_abc import AbstractTransformation
from pymatgen.transformations.standard_transformations import (
    ConventionalCellTransformation,
)
from ase import Atoms


def apply_transformation(tr: AbstractTransformation, atoms: Atoms) -> Atoms:
    adaptor = AseAtomsAdaptor()
    struc = adaptor.get_structure(atoms)
    c = tr.apply_transformation(struc)
    return adaptor.get_atoms(c)


def conventional_cell(atoms: Atoms) -> Atoms:
    return apply_transformation(ConventionalCellTransformation(), atoms)
