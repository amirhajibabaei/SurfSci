# +
from ase.calculators.lammps import convert

from surfsci.lammps.commands import pair_coef_commands, set_commands

topology = {
    "bond": [("O", "H", 1.0)],
    "angle": [("H", "O", "H", 105.26)],
}


# UNITS: real
_lj_params = {
    "O": (0.1553, 3.166),
    "Na": (0.352642, 2.1595),
    "Cl": (0.012785, 4.8305),
    # ramaining is zero
}

# UNITS: real
_charge_params = {"O": -0.8476, "H": 0.4238, "Na": 1, "Cl": -1}


def spce_params(units: str) -> tuple[dict[str, tuple[float, float]], dict[str, float]]:
    """
    returns lj coefs and charges in requested units.

    """
    # lj coefs
    lj_coefs = {}
    for a, (eps, sig) in _lj_params.items():
        eps_ = convert(eps, "energy", "real", units)
        sig_ = convert(sig, "distance", "real", units)
        lj_coefs[a] = (eps_, sig_)

    # charges
    charges = {
        a: convert(q, "charge", "real", units) for a, q in _charge_params.items()
    }
    return lj_coefs, charges


def styles(units: str, cutoff: float = 10.0) -> list[str]:
    """
    cutoff: ASE units (Angstrom)
    """
    cutoff = convert(cutoff, "distance", "ASE", units)

    s = [
        f"pair_style lj/cut/coul/long {cutoff}",
        "pair_modify shift yes mix arithmetic",
        "bond_style  harmonic",
        "angle_style harmonic",
    ]
    return s


def coeffs(units: str, types: dict[str, int], bonds, angles) -> list[str]:
    lj_coefs, charges = spce_params(units)
    oh = bonds[("O", "H")]
    hoh = angles[("H", "O", "H")]
    c = [
        *pair_coef_commands(types, lj_coefs, zero_missing=True),  # type: ignore
        *set_commands("charge", types, charges),
        "kspace_style pppm 0.00001",
        "dielectric 1.0",
        f"bond_coeff {oh} 0 1.000",
        f"angle_coeff {hoh} 0 109.470",
    ]
    return c


def constraints(bonds, angles):
    oh = bonds[("O", "H")]
    hoh = angles[("H", "O", "H")]
    c = [f"fix hoh all rattle 1e-6 100 0 b {oh} a {hoh}"]
    return c
