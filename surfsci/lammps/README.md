<!-- #region -->
# `lammps` subpackage

This subpackage is a wrapper around the original LAMMPS package.
It mainly aims at defining a set of utility commands for automation of workflows.
But we may also define external fixes, etc which can be integrated
with LAMMPS through mechanisms like `fix external`.


## Getting started
A general LAMMPS input script can be executed by
```bash
mpirun -n python -m surfsci.lammps -i in.lammps
```
Under the hood this uses the LAMMPS python interface.
If the input script contains some of extensions which are defined here,
they will be converted into LAMMPS commands.
<!-- #endregion -->
