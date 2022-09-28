# +
"""
This is a work in progress.
"""
from lammps import lammps

import surfsci.lammps.forcefields.demo_nacl as ff
from surfsci.build.nacl import nacl_surface
from surfsci.lammps.atoms import create_atoms


def main() -> None:
    # create NaCl interface with water as Atoms object
    atoms = nacl_surface(indices=(1, 1, 1), size=(6, 6, 3), vacuum=3, h2o=20.0)

    # make LAMMPS object
    cutoff = 10.0
    lmp = lammps()
    lmp.command("units real")
    lmp.command("atom_style full")
    create_atoms(lmp, atoms, box_kw=dict(bond=(1, 2), angle=(1, 1)))
    ff.ff(lmp, atoms, cutoff)
    ff.ff_coef(lmp)

    # write LAMMPS data file
    lmp.command("write_data nacl.dat")

    # Or continue with simulation ...
    # TODO: ...


if __name__ == "__main__":
    main()
