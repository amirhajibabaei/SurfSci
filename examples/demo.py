# +
"""
This is a work in progress.
"""
import surfsci.lammps.forcefields.demo_nacl as ff
from surfsci.build.nacl import nacl_surface
from surfsci.lammps.atoms import atoms_to_lmp


def main() -> None:
    # create NaCl interface with water as Atoms object
    atoms = nacl_surface(indices=(1, 1, 1), size=(6, 6, 3), vacuum=3, h2o=20.0)

    # make LAMMPS object
    cutoff = 10.0
    topology = {
        "bond": [("O", "H", 1.0)],
        "angle": [(0, 0, 105.26)],
    }
    lmp = atoms_to_lmp(atoms, topology=topology)
    ff.ff(lmp, atoms, cutoff)
    ff.ff_coef(lmp)

    # write LAMMPS data file
    lmp.command("write_data nacl.dat")

    # Or continue with simulation ...
    # TODO: ...


if __name__ == "__main__":
    main()
