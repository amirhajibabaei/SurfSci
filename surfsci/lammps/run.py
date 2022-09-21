# +
from typing import Sequence

from lammps import lammps


def run_commands(lmp: lammps, commands: Sequence[str]) -> lammps:
    lmp.commands_list(commands)
    return lmp
