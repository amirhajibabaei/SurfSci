# +
from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.calculators.lammps import convert
from ase.data import atomic_masses, atomic_names
from ase.geometry import wrap_positions
from lammps import lammps

from ..geometry.rotations import PrismRotation
from ..geometry.topology import find_topology


def create_atoms(lmp: lammps, atoms: Atoms, topology: dict = None) -> None:
    """
    prerequisites:
        LAMMPS "units" must be defined before.

    """

    units = lmp.extract_global("units")

    # boundary conditions
    symbol = {True: "p", False: "f"}
    boundary = " ".join([symbol[b] for b in atoms.pbc])
    lmp.command(f"boundary {boundary}")

    # rotate coordinates to LAMMPS "prism" style (most generic)
    # TODO: consider checks for simpler "block" style
    rotation = PrismRotation(atoms.cell)
    cell = rotation(atoms.cell)
    positions = rotation(atoms.positions)
    positions = wrap_positions(positions, cell, atoms.pbc)

    # convert to LAMMPS distance units
    cell = convert(cell, "distance", "ASE", units)
    positions = convert(positions, "distance", "ASE", units)

    # define region
    region_id = "cell"
    lower = cell.flat[[0, 4, 8, 3, 6, 7]]
    prism = "0 {} 0 {} 0 {}  {} {} {}".format(*lower)
    lmp.command(f"region {region_id} prism  {prism} units box")

    # topology
    box_topo = []
    if topology is not None:
        topo = find_topology(atoms, topology)
        for a, b in topology.items():
            c = max([len(c) for c in topo[a]])
            box_topo.append(f"{a}/types {len(b)} extra/{a}/per/atom {c}")
    box_kw = " ".join(box_topo)

    # create box
    unique = np.unique(atoms.get_atomic_numbers())  # auto-sorted
    lmp.command(f"create_box {len(unique)} {region_id} {box_kw}")

    # create atoms
    mapping = {z: i + 1 for i, z in enumerate(unique)}
    types = list(map(mapping.get, atoms.numbers))
    lmp.create_atoms(
        n=len(atoms),
        id=None,
        type=types,
        x=positions.reshape(-1),
        v=None,
        image=None,
        shrinkexceed=False,
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
