import pytest
import sage.all as sa
from kcf import kcf_sage as kcf



def test_L0_L0T():
    A = sa.matrix(sa.SR, [[1]])
    B = sa.matrix(sa.SR, [[0]])
    while True:
        P = sa.random_matrix(sa.ZZ, *A.dimensions()).change_ring(sa.SR)
        if not (P.det().is_zero()):
            break
    A = P.inverse() * A * P
    B = P.inverse() * B * P
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert ((L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero()
            and not L.det().is_zero() and not R.det().is_zero()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_L2_N():
    """
    L_2; N
    """
    A = sa.matrix(sa.SR, [[2, 1, 3], [3, 2, 5], [3, 2, 6]])
    B = sa.matrix(sa.SR, [[1, 1, 2], [1, 1, 2], [1, 1, 3]])
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\nExpected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}"
