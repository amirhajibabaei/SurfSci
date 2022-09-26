# +
from __future__ import annotations

from pathlib import Path

from ase import Atoms
from ase.io import read as read_structure
from mp_api.client import MPRester
from pymatgen.io.ase import AseAtomsAdaptor

# Global:
_root = Path("~/.surfsci/matproj").expanduser()
_api_key = _root.joinpath("api_key.txt")
_cache = _root.joinpath("cache")
_adaptor = AseAtomsAdaptor()


# Special structures
# info -> (mp-id, symmetry, international number, common name)
special_structures: dict[str, list[tuple[str, str, int, str | None]]]
special_structures = {
    "ice": [
        ("mp-32959", "tetragonal", 102, None),
        ("mp-703459", "hexagonal", 185, None),
    ],
    "NaCl": [("mp-22851", "cubic", 221, None), ("mp-22862", "cubic", 225, "rocksalt")],
}


def save_api_key(api_key: str) -> None:
    _root.mkdir(parents=True, exist_ok=True)
    with _api_key.open("w") as of:
        of.write(api_key)


def get_saved_api_key() -> str:
    if not _api_key.exists():
        raise RuntimeError("No saved api-key is found!")
    with _api_key.open("r") as of:
        api_key = of.read()
    return api_key


def get_MPRester(api_key: str | None = None) -> MPRester:
    if api_key is None:
        api_key = get_saved_api_key()
    return MPRester(api_key)


def cache_structure(struc: Atoms, mp_id: str, format: str = "extxyz") -> None:
    folder = _cache.joinpath("struc")
    folder.mkdir(parents=True, exist_ok=True)
    file = folder.joinpath(f"{mp_id}.{format}")
    struc.write(file, format=format)


def cached_structure(mp_id: str, format: str = "extxyz") -> Atoms | None:
    folder = _cache.joinpath("struc")
    file = folder.joinpath(f"{mp_id}.{format}")
    if file.exists():
        return read_structure(file, format=format)
    else:
        return None


def get_structure(
    mp_id: str,
    api_key: str | None = None,
    mp: MPRester | None = None,
    cache: bool = True,
) -> Atoms:
    """
    mp_id: materials project id

    """
    # check cached
    if cache:
        cached = cached_structure(mp_id)
        if cached is not None:
            return cached

    # download
    if mp is None:
        mp = get_MPRester(api_key)
    struc = mp.get_structure_by_material_id(mp_id)
    atoms = _adaptor.get_atoms(struc)

    # cache
    if cache:
        cache_structure(atoms, mp_id)
    return atoms


def find_structure(
    name: str,
    system: str | None = None,
    inum: int | None = None,
    common: str | None = None,
) -> Atoms:
    """
    name: e.g. ice, nacl, etc
    system: crystal system e.g. cubic, hexagonal, etc
    inum: international number (spacegroup)
    common: common name e.g. rocksalt
    """

    def cond(a, b):
        return a == b or b is None

    def match(struc):
        a = cond(struc[1], system)
        b = cond(struc[2], inum)
        c = cond(struc[3], common)
        return a and b and c

    f = list(filter(match, special_structures[name]))
    if len(f) == 0:
        raise RuntimeError("Structure not found!")
    elif len(f) > 1:
        raise RuntimeError("Multiple candidates found!")
    mp_id = f[0][0]
    return get_structure(mp_id)


def _download_special_structures() -> None:
    for system, form in special_structures.items():
        for mp_id, _, _, _ in form:
            get_structure(mp_id)
