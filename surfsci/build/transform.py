# +
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.transformations.standard_transformations import (
    ConventionalCellTransformation,
)
from ase import Atoms

_adaptor = AseAtomsAdaptor()
_conv_tr = ConventionalCellTransformation()


def conventional_cell(atoms: Atoms) -> Atoms:
    struc = _adaptor.get_structure(atoms)
    c = _conv_tr.apply_transformation(struc)
    return _adaptor.get_atoms(c)
