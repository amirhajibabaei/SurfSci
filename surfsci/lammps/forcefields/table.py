# +
from io import StringIO

import numpy as np
from ase.calculators.lammps import convert


def create_table(function, pairs, units_in, units_out=None, file=None, **kwargs):
    if units_out is None:
        units_out = units_in
    table = {}
    for pair, args in pairs.items():
        r, e, f = tabular(function, args, **kwargs)
        r = convert(r, "distance", units_in, units_out)
        e = convert(e, "energy", units_in, units_out)
        f = convert(f, "force", units_in, units_out)
        table[pair] = (r, e, f)
    if file is not None:
        keys = write_table(file, table, units=units_out)
    return table, keys


def tabular(
    function,
    args=None,
    kwargs=None,
    rmin=0.1,
    rmax=10.0,
    dr=0.01,
    cutoff=None,
    shift=True,
):

    if args is None:
        args = ()
    if kwargs is None:
        kwargs = {}

    r = np.arange(rmin, rmax + 1e-8, dr)
    e, f = function(r, *args, **kwargs)

    c = -1
    if cutoff is not None:
        c = np.argmin(abs(r - cutoff))

    if shift:
        e -= e[c]
        f -= f[c]

    if cutoff is not None:
        e[c:] = 0
        f[c:] = 0

    return r, e, f


def write_table(file, tab, units=None):
    with open(file, "w") as of:
        # header
        if units is not None:
            of.write(f"# UNITS: {units}\n")
        # body
        keys = {}
        for pair, (r, e, f) in tab.items():
            N = r.shape[0]
            i = np.arange(1, N + 1)
            data = np.c_[i, r, e, f]
            key = "_".join(pair)
            keys[pair] = (key, N)
            of.write(f"\n{key}")
            of.write(f"\nN {N}\n\n")
            np.savetxt(of, data, fmt=("%10d", "%.18e", "%.18e", "%.18e"))
    return keys


def _read_blocks(path):
    with open(path) as of:
        blocks = [[]]
        for line in of.readlines():
            if line.strip() == "":
                blocks.append([])
            else:
                blocks[-1].append(line)
    blocks = list(filter(lambda b: len(b) > 0, blocks))
    return blocks


def _get_header_sections(blocks):
    nb = len(blocks)
    start = nb % 2
    sections = [(blocks[i], blocks[i + 1]) for i in range(start, nb, 2)]
    if start == 0:
        header = []
    else:
        header = blocks[0]
    return header, sections


def _read_section(section):
    head, data = section
    assert len(head) == 2, "needs generalization!"
    key = head[0].split()[0]
    key = tuple(key.split("_"))
    N, n = head[1].split()
    assert N == "N", "needs generalization!"
    val = np.loadtxt(StringIO("".join(data)))
    assert val.shape[0] == int(n)
    return key, val


def _read_sections(sections):
    out = {}
    for sec in sections:
        key, data = _read_section(sec)
        out[key] = tuple(data[:, k] for k in range(1, 4))
    return out


def _read_table(path):
    blocks = _read_blocks(path)
    header, sections = _get_header_sections(blocks)
    out = _read_sections(sections)
    return out
