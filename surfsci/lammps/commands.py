# +
from __future__ import annotations

import itertools as it
from typing import Any, Sequence


def pair_coef_commands(
    types: dict[str, int],
    coef: dict[str | tuple[str, str], Any],
    explicit_zeros: bool = True,
    pair_style: str | None = None,
):

    if pair_style is None:
        pair_style = ""

    _coef: dict[tuple[int, int], Any] = {}
    null: list[tuple[int, int]] = []
    # same types
    for a in types.keys():
        t = (types[a], types[a])
        if a in coef:
            assert (a, a) not in coef
            _coef[t] = coef[a]
        elif (a, a) in coef:
            assert a not in coef
            _coef[t] = coef[(a, a)]
        else:
            null.append(t)

    # different types
    for (a, b) in it.combinations(types.keys(), 2):
        t = (types[a], types[b])
        if (a, b) in coef:
            assert (b, a) not in coef
            _coef[t] = coef[(a, b)]
        elif (b, a) in coef:
            assert (a, b) not in coef
            _coef[t] = coef[(b, a)]
        else:
            null.append(t)

    # commands
    cmds = []
    for (n, m), c in _coef.items():
        if isinstance(c, Sequence):
            c_seq = c
        else:
            c_seq = (c,)
        c_str = " ".join(map(str, c_seq))
        cmds.append(f"pair_coeff {n} {m} {pair_style} {c_str}")
    if explicit_zeros:
        # c is the last of the previous for-loop
        c_str = " ".join(("0" for _ in c))
        for n, m in null:
            cmds.append(f"pair_coeff {n} {m} {pair_style} {c_str}")
    return cmds


def set_commands(
    quantity: str, types: dict[str, int], values: dict[str, Any]
) -> list[str]:
    return [f"set type {t} {quantity} {values[a]}" for a, t in types.items()]
