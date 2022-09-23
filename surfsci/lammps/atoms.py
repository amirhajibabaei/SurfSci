# +
import numpy as np
from ase import Atoms
from ase.calculators.lammps import convert
from ase.data import atomic_masses, atomic_names
from ase.geometry import wrap_positions
from lammps import lammps


def create_atoms(lmp: lammps, atoms: Atoms) -> None:
    """ """

    units = lmp.extract_global("units")

    # boundary
    symbol = {True: "p", False: "f"}
    boundary = " ".join([symbol[b] for b in atoms.pbc])
    lmp.command(f"boundary {boundary}")

    # define region
    region_id = "region"
    cell = convert(atoms.cell, "distance", "ASE", units)
    if all(atoms.cell.cellpar()[3:] == 90):  # cubic
        region = " ".join(["block"] + [f"{0} {l}" for l in cell.diagonal()])
    else:
        # TODO: generalize
        raise RuntimeError("Non-cubic box is not supported yet!")
    lmp.command(f"region {region_id} {region}")

    # create box
    unique = np.unique(atoms.get_atomic_numbers())  # auto-sorted
    lmp.command(f"create_box {len(unique)} {region_id}")

    # create atoms
    mapping = {z: i + 1 for i, z in enumerate(unique)}
    types = list(map(mapping.get, atoms.numbers))
    _x = convert(atoms.positions, "distance", "ASE", units)
    x = wrap_positions(_x, cell, atoms.pbc).reshape(-1)
    lmp.create_atoms(
        n=len(atoms), id=None, type=types, x=x, v=None, image=None, shrinkexceed=False
    )

    # set masses
    for z, t in mapping.items():
        m = convert(atomic_masses[z], "mass", "ASE", units)
        lmp.command(f"mass {t} {m}")

    # log
    lmp.command("""print "\ntype -> atom" """)
    for z, t in mapping.items():
        lmp.command(f"""print "{t} -> {z} ({atomic_names[z]})" """)

    # set lmp attributes
    lmp._types_mapping = mapping
