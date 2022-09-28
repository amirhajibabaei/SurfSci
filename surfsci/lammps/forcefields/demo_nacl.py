# +
from ase import Atoms
from ase.data import chemical_symbols
from lammps import lammps


def ff(lmp: lammps, atoms: Atoms, cutoff: float) -> None:
    lmp.commands_list(
        [
            f"pair_style  lj/cut/coul/long {cutoff}",
            "pair_modify shift yes mix arithmetic",
            "bond_style harmonic",
            "angle_style harmonic",
        ]
    )
    for a in atoms:
        if a.symbol == "O":
            # LAMMPS ids start from 1
            i = a.index + 1
            lmp.commands_list(
                [
                    f"create_bonds single/bond 1 {i} {i+1}",
                    f"create_bonds single/bond 1 {i} {i+2}",
                    f"create_bonds single/angle 1 {i+1} {i} {i+2}",
                ]
            )


def ff_coef(lmp: lammps) -> None:

    # TODO: generalize to all units
    units = lmp.extract_global("units")
    assert units == "real"

    lj = {
        "O": (0.1553, 3.166),
        "H": (0, 0),
        "Na": (0.352642, 2.1595),
        "Cl": (0.012785, 4.8305),
    }

    charge = {"O": -0.8476, "H": 0.4238, "Na": -1, "Cl": 1}

    for z, t in lmp._types_mapping.items():
        e = chemical_symbols[z]
        eps, sig = lj[e]
        lmp.command(f"pair_coeff {t} {t} {eps} {sig}")
        lmp.command(f"set type {t} charge {charge[e]}")
        lmp.command(f"group {e} type {t}")

    lmp.commands_list(
        [
            "kspace_style pppm 0.00001",
            "dielectric 1.0",
            "bond_coeff 1 100 1.000",
            "angle_coeff 1 300 109.470",
        ]
    )
