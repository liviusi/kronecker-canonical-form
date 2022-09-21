import pytest
import sage.all as sa
from kcf import kcf_sage as kcf


def test_linear_dependence_columns():
    A = sa.matrix(sa.SR, [0, 1, 0])
    B = sa.matrix(sa.SR, [1, 0, 0])
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.is_singular()):
            while True:
                C = sa.random_matrix(sa.ZZ, A.ncols(),
                                     A.ncols()).change_ring(sa.SR)
                if not (C.is_zero()):
                    break
            break
    A = D.inverse() * A * C
    B = D.inverse() * B * C
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert ((L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero()
            and not L.is_singular() and not R.is_singular()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_linear_dependence_rows():
    A = sa.matrix(sa.SR, [0, 1, 0]).transpose()
    B = sa.matrix(sa.SR, [1, 0, 0]).transpose()
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.is_singular()):
            while True:
                C = sa.random_matrix(sa.ZZ, A.ncols(),
                                     A.ncols()).change_ring(sa.SR)
                if not (C.is_zero()):
                    break
            break
    A = D.inverse() * A * C
    B = D.inverse() * B * C
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert ((L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero()
            and not L.is_singular() and not R.is_singular()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_L1_L1T():
    A = sa.matrix(sa.SR, [[0, 1, 0], [0, 0, 0], [0, 0, 1]])
    B = sa.matrix(sa.SR, [[1, 0, 0], [0, 0, 1], [0, 0, 0]])
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.is_singular()):
            while True:
                C = sa.random_matrix(sa.ZZ, A.ncols(),
                                     A.ncols()).change_ring(sa.SR)
                if not (C.is_zero()):
                    break
            break
    A = D.inverse() * A * C
    B = D.inverse() * B * C
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert ((L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero()
            and not L.is_singular() and not R.is_singular()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_L2():
    A = sa.matrix(sa.SR, [[0, 1], [0, 0]])
    B = sa.matrix(sa.SR, [[1, 0], [0, 1]])
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.is_singular()):
            while True:
                C = sa.random_matrix(sa.ZZ, A.ncols(),
                                     A.ncols()).change_ring(sa.SR)
                if not (C.is_zero()):
                    break
            break
    A = D.inverse() * A * C
    B = D.inverse() * B * C
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert ((L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero()
            and not L.is_singular() and not R.is_singular()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_L2T():
    A = sa.matrix(sa.SR, [[0, 1], [0, 0]]).transpose()
    B = sa.matrix(sa.SR, [[1, 0], [0, 1]]).transpose()
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.is_singular()):
            while True:
                C = sa.random_matrix(sa.ZZ, A.ncols(),
                                     A.ncols()).change_ring(sa.SR)
                if not (C.is_zero()):
                    break
            break
    A = D.inverse() * A * C
    B = D.inverse() * B * C
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert ((L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero()
            and not L.is_singular() and not R.is_singular()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_L2_L2T():
    A = sa.matrix(sa.SR, [0, 1, 0, 0])
    A = sa.block_diagonal_matrix([A, A.transpose()])
    B = sa.matrix(sa.SR, [1, 0, 0, 0])
    B = sa.block_diagonal_matrix([B, B.transpose()])
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.is_singular()):
            while True:
                C = sa.random_matrix(sa.ZZ, A.ncols(),
                                     A.ncols()).change_ring(sa.SR)
                if not (C.is_zero()):
                    break
            break
    A = D.inverse() * A * C
    B = D.inverse() * B * C
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert ((L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero()
            and not L.is_singular() and not R.is_singular()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_L3():
    A = sa.matrix(sa.SR, [[0, 1, 0], [0, 0, 1]])
    B = sa.matrix(sa.SR, [[1, 0, 0], [0, 1, 0]])
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.is_singular()):
            while True:
                C = sa.random_matrix(sa.ZZ, A.ncols(),
                                     A.ncols()).change_ring(sa.SR)
                if not (C.is_zero()):
                    break
            break
    A = D.inverse() * A * C
    B = D.inverse() * B * C
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert ((L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero()
            and not L.is_singular() and not R.is_singular()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_L4_two_rows():
    A = sa.matrix(sa.SR, [[0, 1, 0, 0], [0, 0, 1, 0]])
    B = sa.matrix(sa.SR, [[1, 0, 0, 0], [0, 1, 0, 0]])
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.is_singular()):
            while True:
                C = sa.random_matrix(sa.ZZ, A.ncols(),
                                     A.ncols()).change_ring(sa.SR)
                if not (C.is_zero()):
                    break
            break
    A = D.inverse() * A * C
    B = D.inverse() * B * C
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert ((L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero()
            and not L.is_singular() and not R.is_singular()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_L4T_two_rows():
    A = sa.matrix(sa.SR, [[0, 1, 0, 0], [0, 0, 1, 0]]).transpose()
    B = sa.matrix(sa.SR, [[1, 0, 0, 0], [0, 1, 0, 0]]).transpose()
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.is_singular()):
            while True:
                C = sa.random_matrix(sa.ZZ, A.ncols(),
                                     A.ncols()).change_ring(sa.SR)
                if not (C.is_zero()):
                    break
            break
    A = D.inverse() * A * C
    B = D.inverse() * B * C
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert ((L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero()
            and not L.is_singular() and not R.is_singular()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_L4():
    A = sa.matrix(sa.SR, [[0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]])
    B = sa.matrix(sa.SR, [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0]])
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.is_singular()):
            while True:
                C = sa.random_matrix(sa.ZZ, A.ncols(),
                                     A.ncols()).change_ring(sa.SR)
                if not (C.is_zero()):
                    break
            break
    A = D.inverse() * A * C
    B = D.inverse() * B * C
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert ((L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero()
            and not L.is_singular() and not R.is_singular()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_L4T():
    A = sa.matrix(sa.SR, [[0, 1, 0, 0, 0], [0, 0, 1, 0, 0],
                          [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]]).transpose()
    B = sa.matrix(sa.SR, [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0],
                          [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]]).transpose()
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.is_singular()):
            while True:
                C = sa.random_matrix(sa.ZZ, A.ncols(),
                                     A.ncols()).change_ring(sa.SR)
                if not (C.is_zero()):
                    break
            break
    A = D.inverse() * A * C
    B = D.inverse() * B * C
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert ((L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero()
            and not L.is_singular() and not R.is_singular()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")
