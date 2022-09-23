# +
"""
Row convention:
    The x, y, z coordinates of a vector in 3d
    are assumed to be stored in rows of a numpy
    arrays with the following layout:
    R = [
            [x1, y1, z1],
            [x2, y2, z2],
            ...
        ]
    This is chosen for consistency with other
    related python packages such as ase, etc.
    But one has to keep in mind that in algebraic
    formulas the vector are usually assumed as
    columns. Therefore the code will look slightly
    different from algebraic formulas. For instance
    rotation of coordinates R by a matrix M becomes:
        rotated-R = (M @ R.T).T = R @ M.T

"""
import abc
from math import pi

import numpy as np
from numpy import cross, dot
from numpy.linalg import det, norm


class Rotation:
    """
    Callables which rotate a set of 3d vectors.
    """

    @abc.abstractmethod
    def __call__(self, xyz: np.ndarray) -> np.ndarray:
        ...


class AxisRotation(Rotation):
    """
    Rotations in 3d Cartesian coordinates defined by
    a rotation "axis" and an "angle".
    """

    def __init__(self, axis: tuple[float, float, float], angle: float) -> None:
        """
        axis: a rotation axis
        angle: rotaion angle (Radian)
        """
        n = np.asarray(axis)
        n = n / norm(n)
        a = np.cos(angle / 2)
        b, c, d = -n * np.sin(angle / 2)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        self._rot_t = np.transpose(
            np.array(
                [
                    [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                    [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                    [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
                ]
            )
        )

    def __call__(self, xyz: np.ndarray) -> np.ndarray:
        return xyz @ self._rot_t


class PrismRotation(Rotation):
    """
    Defines a rotation which maps 3 right-handed
    generic base vectors A, B, C into new ones
    a, b, c where
        a = [ax,  0,  0]
        b = [bx, by,  0]
        c = [cx, cy, cz]

    For more information, see:
        https://docs.lammps.org/Howto_triclinic.html

    """

    def __init__(self, basis: np.ndarray) -> None:
        """
        basis: a 3x3 array where each row is a base
                vector for a right-handed system
        """
        # TODO: raise error if not right-handed
        a, b, c = basis
        a_ = a / norm(a)
        ab = cross(a, b)
        ab_ = ab / norm(ab)
        prism = np.zeros((3, 3))
        prism.flat[[0, 4, 8, 1, 2, 5]] = [
            norm(a),  # ax -> 0
            norm(cross(a_, b)),  # by -> 4
            dot(ab_, c),  # cz -> 8
            dot(a_, b),  # bx -> 1
            dot(a_, c),  # cx -> 2
            dot(cross(ab_, a_), c),  # cy -> 5
        ]
        recip = np.array([cross(b, c), cross(c, a), cross(a, b)]) / det(prism)
        self._rot_t = (prism @ recip).T

    def __call__(self, xyz: np.ndarray) -> np.ndarray:
        return xyz @ self._rot_t


def test_AxisRotation() -> bool:
    axis = (0.0, 0.0, 1.0)
    angle = pi / 2
    rot = AxisRotation(axis, angle)
    xyz_in = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0],
        ]
    )
    xyz_out = np.array(
        [
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, -1.0, 0.0],
        ]
    )
    # one at a time
    for a, b in zip(xyz_in, xyz_out):
        assert np.allclose(rot(a), b)
    # collective
    assert np.allclose(rot(xyz_in), xyz_out)
    return True


def test_PrismRotation() -> bool:
    basis = np.random.uniform(size=(3, 3))
    new = PrismRotation(basis)(basis)
    assert np.allclose(new.flat[[1, 2, 5]], 0)
    assert np.allclose(basis @ basis.T, new @ new.T)
    return True


if __name__ == "__main__":
    test_AxisRotation()
    test_PrismRotation()
