# +
import argparse

from lammps import lammps

from ._io import get_input_commands
from .run import run_commands


def main() -> None:
    # i. command line arguments
    parser = argparse.ArgumentParser(description="LAMMPS runner")
    parser.add_argument("-i", "--input", default=None, type=str, help="input file")
    args = parser.parse_args()
    if args.input is None:
        # TODO: make it consistent with LAMMPS
        raise RuntimeError("An input scrip is needed!")

    # ii. parse and run the input script
    commands = get_input_commands(args.input)
    lmp = lammps()
    lmp = run_commands(lmp, commands)


if __name__ == "__main__":
    main()
