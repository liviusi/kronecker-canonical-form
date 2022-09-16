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


def test_L1_L1T_N2():
    A = sa.matrix(sa.SR, [[0, 1]])
    A = sa.block_diagonal_matrix([A, A.transpose(),
                                  sa.identity_matrix(2)])
    B = sa.matrix(sa.SR, [[1, 0]])
    B = sa.block_diagonal_matrix([B, B.transpose(),
                                  sa.matrix(sa.SR, [[0, 1], [0, 0]])])
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


def test_L1_L1_N():
    A = sa.matrix(sa.SR, [[0, 1]])
    A = sa.block_diagonal_matrix([A, A, sa.identity_matrix(1)])
    B = sa.matrix(sa.SR, [[1, 0]])
    B = sa.block_diagonal_matrix([B, B, sa.matrix(sa.SR, [0])])
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.det().is_zero()):
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
            and not L.det().is_zero() and not R.det().is_zero()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_L1_L2T_N():
    L2_A = sa.matrix(sa.SR, [[0, 1, 0], [0, 0, 1]])
    L2_B = sa.matrix(sa.SR, [[1, 0, 0], [0, 1, 0]])
    A = sa.matrix(sa.SR, [[0, 1]])
    A = sa.block_diagonal_matrix([A, L2_A.transpose(),
                                  sa.identity_matrix(1)])
    B = sa.matrix(sa.SR, [[1, 0]])
    B = sa.block_diagonal_matrix([B, L2_B.transpose(),
                                  sa.matrix(sa.SR, [0])])
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.det().is_zero()):
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
            and not L.det().is_zero() and not R.det().is_zero()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_L2_L1T_J():
    L2_A = sa.matrix(sa.SR, [[0, 1, 0], [0, 0, 1]])
    L2_B = sa.matrix(sa.SR, [[1, 0, 0], [0, 1, 0]])
    A = sa.block_diagonal_matrix([L2_A,
                                  sa.matrix(sa.SR, [[0, 1]]).transpose(),
                                  sa.matrix(sa.SR, [42])])
    B = sa.block_diagonal_matrix([L2_B,
                                  sa.matrix(sa.SR, [[1, 0]]).transpose(),
                                  sa.identity_matrix(1)])
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.det().is_zero()):
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
            and not L.det().is_zero() and not R.det().is_zero()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_L3_N2_J2():
    L3_A = sa.matrix(sa.SR, [[0, 1, 0, 0], [0, 0, 1, 0]])
    L3_B = sa.matrix(sa.SR, [[1, 0, 0, 0], [0, 1, 0, 0]])
    A = sa.block_diagonal_matrix([L3_A,
                                  sa.identity_matrix(2),
                                  sa.matrix(sa.SR, [[42, 1], [0, 42]])])
    B = sa.block_diagonal_matrix([L3_B,
                                  sa.matrix(sa.SR, [[0, 1], [0, 0]]),
                                  sa.identity_matrix(2)])
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.det().is_zero()):
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
            and not L.det().is_zero() and not R.det().is_zero()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_L3T_N2_J2():
    L3_A = sa.matrix(sa.SR, [[0, 1, 0, 0], [0, 0, 1, 0]])
    L3_B = sa.matrix(sa.SR, [[1, 0, 0, 0], [0, 1, 0, 0]])
    A = sa.block_diagonal_matrix([L3_A,
                                  sa.identity_matrix(2),
                                  sa.matrix(sa.SR, [[42, 1],
                                                    [0, 42]])]).transpose()
    B = sa.block_diagonal_matrix([L3_B,
                                  sa.matrix(sa.SR, [[0, 1], [0, 0]]),
                                  sa.identity_matrix(2)]).transpose()
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.det().is_zero()):
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
            and not L.det().is_zero() and not R.det().is_zero()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_L3_N3_J2_F3():
    L3_A = sa.matrix(sa.SR, [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    L3_B = sa.matrix(sa.SR, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
    A = sa.block_diagonal_matrix([L3_A,
                                  sa.identity_matrix(2),
                                  sa.matrix(sa.SR, [[42, 1],
                                                    [0, 42]])]).transpose()
    B = sa.block_diagonal_matrix([L3_B,
                                  sa.matrix(sa.SR, [[0, 1], [0, 0]]),
                                  sa.identity_matrix(2)]).transpose()
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.det().is_zero()):
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
            and not L.det().is_zero() and not R.det().is_zero()), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")