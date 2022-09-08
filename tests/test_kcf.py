import pytest
import sage.all as sa
from kcf import kcf_sage as kcf


def test_one():
    """
    L_1; L_1^T
    """
    A = sa.matrix(sa.SR, [[0, 1, 0], [0, 0, 0], [0, 0, 1]])
    B = sa.matrix(sa.SR, [[1, 0, 0], [0, 0, 1], [0, 0, 0]])
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\nExpected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}"

def test_two():
    """
    L_0; L_0^T
    """
    A = sa.matrix(sa.SR, [[0]])
    B = sa.matrix(sa.SR, [[1]])
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\nExpected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}"

def test_three():
    """
    """
    A = sa.matrix(sa.SR, [[1]])
    B = sa.matrix(sa.SR, [[0]])
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\nExpected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}"

def test_four():
    """
    L_2; N
    """
    A = sa.matrix(sa.SR, [[2, 1, 3], [3, 2, 5], [3, 2, 6]])
    B = sa.matrix(sa.SR, [[1, 1, 2], [1, 1, 2], [1, 1, 3]])
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\nExpected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}"

def test_five():
    """
    """
    A = sa.matrix(sa.SR, [[0, 3, -1, 4], [2, -108, 0, 1], [-2, -1, -1, 1], [4, 1, -4, 1]])
    B = sa.matrix(sa.SR, [[0, 0, 0, -1], [0, 2, 0, 0], [0, -2, 0, -1], [0, 4, 0, -4]])
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\nExpected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}"