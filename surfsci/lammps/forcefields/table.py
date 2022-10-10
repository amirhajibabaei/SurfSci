# +
from io import StringIO

import numpy as np


def write_table(file, tab, units=None):
    with open(file, "w") as of:
        # header
        if units is not None:
            of.write(f"# UNITS: {units}\n")
        # body
        for pair, (r, e, f) in tab.items():
            N = r.shape[0]
            i = np.arange(1, N + 1)
            data = np.c_[i, r, e, f]
            of.write(f"\n{pair}")
            of.write(f"\nN {N}\n\n")
            np.savetxt(of, data, fmt=("%10d", "%.18e", "%.18e", "%.18e"))


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
    N, n = head[1].split()
    assert N == "N", "needs generalization!"
    val = np.loadtxt(StringIO("".join(data)))
    assert val.shape[0] == int(n)
    return key, val


def _read_sections(sections):
    out = {}
    for sec in sections:
        key, data = _read_section(sec)
        out[key] = data
    return out


def _read_table(path):
    blocks = _read_blocks(path)
    header, sections = _get_header_sections(blocks)
    out = _read_sections(sections)
    return out
