# +
import numpy as np


def in_range(f: float, a: float, b: float) -> bool:
    if f >= a and f <= b:
        return True
    else:
        return False


def in_range_rtol(f: float, g: float, rtol: float) -> bool:
    return in_range(f, g * (1 - rtol), g * (1 + rtol))


def array_is_tril(arr: np.ndarray) -> bool:
    return np.allclose(arr[np.triu_indices_from(arr, k=1)], 0.0)


def array_is_triu(arr: np.ndarray) -> bool:
    return np.allclose(arr[np.tril_indices_from(arr, k=-1)], 0.0)


def array_is_diag(arr: np.ndarray) -> bool:
    return array_is_triu(arr) and array_is_tril(arr)
