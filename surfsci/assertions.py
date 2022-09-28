# +
import numpy as np


def array_is_tril(arr: np.ndarray) -> bool:
    return np.allclose(arr[np.triu_indices_from(arr, k=1)], 0.0)


def array_is_triu(arr: np.ndarray) -> bool:
    return np.allclose(arr[np.tril_indices_from(arr, k=-1)], 0.0)


def array_is_diag(arr: np.ndarray) -> bool:
    return array_is_triu(arr) and array_is_tril(arr)
