import pytest
import sage.all as sa
from kcf import kcf_sage as kcf


def test_N():
    A = sa.identity_matrix(1).change_ring(sa.SR)
    B = sa.matrix(sa.SR, [0])
    while True:
        P = sa.random_matrix(sa.ZZ, *A.dimensions()).change_ring(sa.SR)
        if not (P.det().is_zero()):
            break
    A = P.inverse() * A * P
    B = P.inverse() * B * P
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_J():
    A = sa.matrix(sa.SR, [42])
    B = sa.matrix(sa.SR, [1])
    while True:
        P = sa.random_matrix(sa.ZZ, *A.dimensions()).change_ring(sa.SR)
        if not (P.det().is_zero()):
            break
    A = P.inverse() * A * P
    B = P.inverse() * B * P
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_N_N():
    A = sa.identity_matrix(2).change_ring(sa.SR)
    B = sa.matrix(sa.SR, [[0, 1], [0, 0]])
    while True:
        P = sa.random_matrix(sa.ZZ, *A.dimensions()).change_ring(sa.SR)
        if not (P.det().is_zero()):
            break
    A = P.inverse() * A * P
    B = P.inverse() * B * P
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_N_J():
    A = sa.matrix(sa.SR, [[1, 0], [0, 42]])
    B = sa.matrix(sa.SR, [[0, 1], [0, 1]])
    while True:
        P = sa.random_matrix(sa.ZZ, *A.dimensions()).change_ring(sa.SR)
        if not (P.det().is_zero()):
            break
    A = P.inverse() * A * P
    B = P.inverse() * B * P
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_N_N_N():
    A = sa.identity_matrix(3).change_ring(sa.SR)
    B = sa.matrix(sa.SR, [[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    while True:
        P = sa.random_matrix(sa.ZZ, *A.dimensions()).change_ring(sa.SR)
        if not (P.det().is_zero()):
            break
    A = P.inverse() * A * P
    B = P.inverse() * B * P
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_N_N_J():
    A = sa.matrix(sa.SR, [[1, 0, 0], [0, 1, 0], [0, 0, 42]])
    B = sa.matrix(sa.SR, [[0, 1, 0], [0, 0, 1], [0, 0, 1]])
    while True:
        P = sa.random_matrix(sa.ZZ, *A.dimensions()).change_ring(sa.SR)
        if not (P.det().is_zero()):
            break
    A = P.inverse() * A * P
    B = P.inverse() * B * P
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_N_J_J():
    A = sa.matrix(sa.SR, [[1, 0, 0], [0, 42, 1], [0, 0, 42]])
    B = sa.matrix(sa.SR, [[0, 1, 0], [0, 1, 0], [0, 0, 1]])
    while True:
        P = sa.random_matrix(sa.ZZ, *A.dimensions()).change_ring(sa.SR)
        if not (P.det().is_zero()):
            break
    A = P.inverse() * A * P
    B = P.inverse() * B * P
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_J_J_J():
    A = sa.matrix(sa.SR, [[42, 0, 0], [0, 5, 0], [0, 0, 2]])
    B = sa.matrix(sa.SR, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    while True:
        P = sa.random_matrix(sa.ZZ, *A.dimensions()).change_ring(sa.SR)
        if not (P.det().is_zero()):
            break
    A = P.inverse() * A * P
    B = P.inverse() * B * P
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_N_N_N_N():
    A = sa.identity_matrix(4).change_ring(sa.SR)
    B = sa.matrix(sa.SR, [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, 0]])
    while True:
        P = sa.random_matrix(sa.ZZ, *A.dimensions()).change_ring(sa.SR)
        if not (P.det().is_zero()):
            break
    A = P.inverse() * A * P
    B = P.inverse() * B * P
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_N_N_N_J():
    A = sa.matrix(sa.SR, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 42]])
    B = sa.matrix(sa.SR, [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, 1]])
    while True:
        P = sa.random_matrix(sa.ZZ, *A.dimensions()).change_ring(sa.SR)
        if not (P.det().is_zero()):
            break
    A = P.inverse() * A * P
    B = P.inverse() * B * P
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_N_N_J_J():
    A = sa.matrix(sa.SR, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 42, 1], [0, 0, 0, 42]])
    B = sa.matrix(sa.SR, [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    while True:
        P = sa.random_matrix(sa.ZZ, *A.dimensions()).change_ring(sa.SR)
        if not (P.det().is_zero()):
            break
    A = P.inverse() * A * P
    B = P.inverse() * B * P
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_N_J_J_J():
    A = sa.matrix(sa.SR, [[1, 0, 0, 0], [0, 42, 1, 0], [0, 0, 42, 1], [0, 0, 0, 42]])
    B = sa.matrix(sa.SR, [[0, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    while True:
        P = sa.random_matrix(sa.ZZ, *A.dimensions()).change_ring(sa.SR)
        if not (P.det().is_zero()):
            break
    A = P.inverse() * A * P
    B = P.inverse() * B * P
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")


def test_J_J_J_J():
    A = sa.matrix(sa.SR, [[42, 1, 0, 0], [0, 42, 1, 0], [0, 0, 42, 1], [0, 0, 0, 42]])
    B = sa.matrix(sa.SR, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    while True:
        P = sa.random_matrix(sa.ZZ, *A.dimensions()).change_ring(sa.SR)
        if not (P.det().is_zero()):
            break
    A = P.inverse() * A * P
    B = P.inverse() * B * P
    (L, R), (KCF_A, KCF_B) = kcf.kronecker_canonical_form(A, B, True)
    assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero(), (
        f"Gotten:\n{kcf.stringify_pencil(KCF_A, KCF_B)}\n",
        f"Expected:\n{kcf.stringify_pencil(L*A*R, L*B*R)}")
