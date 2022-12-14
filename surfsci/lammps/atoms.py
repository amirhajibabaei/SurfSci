# +
from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.calculators.lammps import convert
from ase.data import chemical_symbols
from ase.geometry import wrap_positions
from lammps import lammps

from ..geometry.rotations import PrismRotation
from ..geometry.topology import find_topology
from .commands import mass_commands


class Wlammps:
    def __init__(self, lmp: lammps) -> None:
        self.lmp = lmp
        self.create_atoms = lmp.create_atoms
        self.extract_global = lmp.extract_global

        self.history: list[str | list[str]] = []
        self._num_types: dict[int, int]
        self._chem_types: dict[str, int]
        self._topology: dict | None

    def commands_list(self, cmds: list[str], save=True) -> None:
        if save:
            self.history.append(cmds)
        self.lmp.commands_list(cmds)

    def command(self, cmd: str, save=True) -> None:
        if save:
            self.history.append([cmd])
        self.lmp.command(cmd)

    def forcefield(self, styles, coeffs, consts):
        units = self.extract_global("units")
        self.commands_list(styles(units))
        self.command("# write_data")
        bonds = {b[:2]: self._type(*b[:2]) for b in self._topology["bond"]}
        angles = {b[:3]: self._type(*b[:3]) for b in self._topology["angle"]}
        self.commands_list(coeffs(units, self._chem_types, bonds, angles))
        self.commands_list(consts(bonds, angles))

    def generate_input_files(self):
        with open("in.lammps", "w") as of:
            for b in self.history:
                of.write("\n")
                for l in b:
                    if l.startswith("# write_data"):
                        self.lmp.command("write_data box.in")
                        of.write("\nread_data box.in\n")
                    else:
                        of.write(l + "\n")

    def _bond_type(self, a: str, b: str) -> int:
        if self._topology is not None:
            for i, x in enumerate(self._topology["bond"]):
                if set(x[:-1]) == set((a, b)):
                    return i + 1
        return 0

    def _angle_type(self, a: str, b: str, c: str) -> int:
        if self._topology is not None:
            for i, x in enumerate(self._topology["angle"]):
                if tuple(x[:-1]) == (a, b, c):
                    return i + 1
        return 0

    def _type(self, *at: str) -> int:
        if len(at) == 1:
            return self._chem_types[at[0]]
        elif len(at) == 2:
            return self._bond_type(*at)
        elif len(at) == 3:
            return self._angle_type(*at)
        else:
            return 0


def atoms_to_lmp(atoms: Atoms, topology: dict = None, units="real") -> Wlammps:
    """
    TODO: doc

    """

    # units = lmp.extract_global("units")
    lmp = Wlammps(lammps())

    # units and style
    lmp.commands_list(
        [
            f"units {units}",
            "atom_style full",
        ]
    )

    # boundary conditions
    symbol = {True: "p", False: "f"}
    boundary = " ".join([symbol[b] for b in atoms.pbc])
    lmp.command(f"boundary {boundary}", save=False)

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
    lmp.command(f"region {region_id} prism  {prism} units box", save=False)

    # topology
    topologies = {}
    _box_extra = []
    if topology is not None:
        topologies = find_topology(atoms, topology)
        for a, b in topology.items():
            c = max([len(at) for at in topologies[a]])
            _box_extra.append(f"{a}/types {len(b)} extra/{a}/per/atom {c}")
    box_extra = " ".join(_box_extra)

    # create box
    unique = np.unique(atoms.get_atomic_numbers())  # auto-sorted
    lmp.command(f"create_box {len(unique)} {region_id} {box_extra}", save=False)

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

    # create topologies
    for style, atoms_ in topologies.items():
        cmd = []
        for atom_ in atoms_:
            for top in atom_:
                type_ = top[0] + 1
                ind = " ".join(map(str, [i + 1 for i in top[1:]]))  # TODO: map indices
                cmd.append(f"create_bonds single/{style} {type_} {ind} special no")
        if len(cmd) > 0:
            lmp.commands_list(cmd[:-1], save=False)
            lmp.command(cmd[-1].replace("special no", "special yes"), save=False)

    # set lmp attributes
    lmp._num_types = mapping
    lmp._chem_types = {chemical_symbols[z]: t for z, t in mapping.items()}
    lmp._topology = topology

    # set masses
    lmp.commands_list(mass_commands(lmp._chem_types, units), save=False)
    lmp.commands_list([f"# {e}: {t}" for e, t in lmp._chem_types.items()])

    return lmp
