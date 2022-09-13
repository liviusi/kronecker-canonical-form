import sage.all as sa
from random import randint


def _complete_to_a_basis(M: sa.sage.matrix) -> sa.sage.matrix:
    """
    Completes the given matrix to a basis by augmenting
    it with the identity matrix and recovering the pivots.
    """
    identity_mat = sa.matrix(M.base_ring(), [
        [1 if i == j else 0 for i in range(M.nrows())]
        for j in range(M.nrows())])
    new_M = M.augment(identity_mat)
    return new_M[:, new_M.pivots()]


def _build_coefficient_matrix(
        A: sa.sage.matrix,
        B: sa.sage.matrix,
        epsilon: int) -> sa.sage.matrix:
    """
    Given two matrices (A, B) of the same size returns a matrix of
    the following form:
    [A   0     0]
    [B   A      ]
    [0   B      ]
    [          A]
    [0 0  ...  B]
    of size (epsilon + (2 x epsilon + 1)).
    """
    assert A.ncols() == B.ncols() and A.nrows() == B.nrows()
    assert A.base_ring() == B.base_ring()
    if epsilon <= 0:
        return sa.matrix(A.base_ring(), [])
    else:
        rows = []
        if epsilon == 1:
            rows.append([A])
            rows.append([B])
        else:
            for i in range(epsilon):
                rows.append([A if i == j else
                             (B if i == j + 1 else 0)
                             for j in range(epsilon)])
            rows.append([B if j == epsilon - 1 else 0 for j in range(epsilon)])
        return sa.block_matrix(A.base_ring(), rows)


def _compute_lowest_degree_polynomial(
        A: sa.sage.matrix,
        B: sa.sage.matrix) -> tuple[int, list[sa.vector]]:
    """
    Returns the pair (degree, v) with v representing the
    polynomial of minimum degree in the kernel of the pencil
    (A + tB)x = 0
    and degree its degree.
    """
    assert A.nrows() == B.nrows() and B.ncols() == A.ncols()
    assert A.base_ring() == B.base_ring()
    degree = 0
    coefficient_matrix_kernel = None
    while True:
        coefficient_matrix = _build_coefficient_matrix(A, B, degree + 1)
        coefficient_matrix_kernel = (
            coefficient_matrix.right_kernel().basis_matrix()
        )
        if coefficient_matrix_kernel:
            break
        degree += 1
    coefficient_matrix_kernel = coefficient_matrix_kernel.transpose()
    indices = [i for i in range(0,
                                coefficient_matrix_kernel.nrows(),
                                coefficient_matrix_kernel.nrows() // (
                                    degree + 1))]
    coefficient_matrix_kernel.subdivide(indices, None)
    result = []
    for index in indices:
        col = []
        for i in range(index,
                       index + coefficient_matrix_kernel.nrows() // (
                           degree + 1)):
            col.append(coefficient_matrix_kernel[i, 0])
        result.append(sa.vector(A.base_ring(), col))

    return (degree, result)


def _reduction_theorem(
        A: sa.sage.matrix,
        B: sa.sage.matrix,
        transformation=False):
    """
    Reduces a pencil (A + tB)x = 0 to a canonical form
    A_tilde = P^-1 * A * Q =
    [0 1 0 ... 0 |        ]
    [0 0 1 ... 0 |   D    ]
    [0 0 0 ... 1 |        ]
    [------------|------- ]
    [     0      | A_STAR ]

    B_tilde = P^-1 * B * Q =
    [1 0 0 ... 0 |        ]
    [0 1 0 ... 0 |   F    ]
    [0 0 ... 1 0 |        ]
    [------------|------- ]
    [     0      | B_STAR ]

    and returns A_tilde, B_tilde.

    If transformation is toggled on,
    it returns the tuple (P, Q), (A_tilde, B_tilde) with
    P, Q such that
    P*A*Q = A_tilde,
    P*B*Q = B_tilde.
    """
    assert A.nrows() == B.nrows() and A.ncols() == B.ncols()
    EMPTY_MATRIX = sa.matrix(A.base_ring(), [])
    degree, polynomial = _compute_lowest_degree_polynomial(A, -B)

    V = sa.matrix(polynomial[0].base_ring(), polynomial).transpose()
    Q = _complete_to_a_basis(V)
    P = _complete_to_a_basis(A * V)
    A_tilde = P**-1 * A * Q
    B_tilde = P**-1 * B * Q

    if degree == 0:
        if transformation:
            return ((P**-1, Q),
                    (sa.matrix(sa.SR, [
                        0 for _ in range(A_tilde.nrows())]).transpose(),
                     A_tilde.submatrix(0, 1)),
                    (sa.matrix(sa.SR, [
                        0 for _ in range(A_tilde.nrows())]).transpose(),
                     B_tilde.submatrix(0, 1)))
        else:
            return ((sa.matrix(sa.SR, [
                0 for _ in range(A_tilde.nrows())]).transpose(),
                     A_tilde.submatrix(0, 1)),
                    (sa.matrix(sa.SR, [
                        0 for _ in range(A_tilde.nrows())]).transpose(),
                     B_tilde.submatrix(0, 1)))

    L_A, L_B = [], []
    # now partitioning A_tilde, B_tilde
    for i in range(degree):
        row_L_A = []
        row_L_B = []
        for j in range(degree + 1):
            row_L_A.append(A_tilde[i, j])
            row_L_B.append(B_tilde[i, j])
        L_A.append(row_L_A)
        L_B.append(row_L_B)
    L_A = sa.matrix(A_tilde.base_ring(), L_A)
    L_B = sa.matrix(B_tilde.base_ring(), L_B)

    # (A + tB)x = 0 is now in KCF.
    if A.ncols() == L_A.ncols() and A.nrows() == L_B.nrows():
        if not transformation:
            return ((L_A, EMPTY_MATRIX),
                    (L_B, EMPTY_MATRIX))
        else:
            return ((P**-1, Q), (L_A, EMPTY_MATRIX), (L_B, EMPTY_MATRIX))

    D = []
    for i in range(L_A.nrows()):
        row = []
        for j in range(L_A.ncols(), A.ncols()):
            row.append(A_tilde[i, j])
        D.append(row)
    D = sa.matrix(A.base_ring(), D)

    A_STAR = []
    for i in range(L_A.nrows(), A.nrows()):
        row = []
        for j in range(L_A.ncols(), A.ncols()):
            row.append(A_tilde[i, j])
        A_STAR.append(row)
    A_STAR = sa.matrix(A.base_ring(), A_STAR)

    F = []
    for i in range(L_B.nrows()):
        row = []
        for j in range(L_B.ncols(), B.ncols()):
            row.append(B_tilde[i, j])
        F.append(row)
    F = sa.matrix(B.base_ring(), F)

    B_STAR = []
    for i in range(L_B.nrows(), B.nrows()):
        row = []
        for j in range(L_B.ncols(), B.ncols()):
            row.append(B_tilde[i, j])
        B_STAR.append(row)
    B_STAR = sa.matrix(B.base_ring(), B_STAR)

    W = Z = Y = None

    if degree - 1 > 0:
        M = _build_coefficient_matrix(A_STAR, B_STAR, degree-1)
        if degree > 1:
            sign = -1
            W = []
            for i in range(F.nrows() - 1):
                for j in range(F.ncols()):
                    if j % (A.ncols() - degree - 1) == 0:
                        sign *= -1
                    W.append(sign * (F[i+1, j] - D[i, j]))
            W = sa.matrix(A.base_ring(), [W])
            Z = M.solve_left(W)
            assert Z.ncols() % degree == 0
            Y = []
            row = []
            sign = 1
            for i in range(Z.ncols()):
                if i != 0 and i % (Z.ncols() // degree) == 0:
                    sign *= -1
                    Y.append(row)
                    row = []
                row.append(sign * Z[0, i])
            if row:
                Y.append(row)
            Y = sa.matrix(A.base_ring(), Y)
    else:
        Y = sa.matrix(A.base_ring(), [[1 for _ in range(B_STAR.nrows())]])

    X = [[None for _ in range(A_STAR.ncols())] for _ in range(degree + 1)]
    for i in range(degree):
        for j in range(A_STAR.ncols()):
            X[i][j] = F[i, j] + Y.row(i) * B_STAR.column(j)
    for j in range(A_STAR.ncols()):
        X[degree][j] = D[degree - 1, j] + Y.row(degree - 1) * A_STAR.column(j)

    X = sa.matrix(A.base_ring(), X)

    assert ((D + Y * A_STAR - L_A * X).is_zero()) and (
        (F + Y * B_STAR - L_B * X).is_zero())

    if not transformation:
        return ((L_A, A_STAR), (L_B, B_STAR))
    else:
        left_M = sa.block_matrix(sa.SR, [[sa.identity_matrix(degree), Y],
                                         [0, sa.identity_matrix(
                                            A_tilde.nrows() - degree)]])
        right_M = sa.block_matrix(sa.SR, [[sa.identity_matrix(degree + 1), -X],
                                          [0, sa.identity_matrix(
                                            A_tilde.ncols() - degree - 1)]])
        return ((left_M * P**-1, Q * right_M), (L_A, A_STAR), (L_B, B_STAR))


def _reduce_regular_pencil(A: sa.sage.matrix,
                           B: sa.sage.matrix,
                           transformation=False) -> tuple[tuple[
                                sa.sage.matrix, sa.sage.matrix],
                                                          tuple[sa.sage.matrix,
                                                                sa.sage.matrix]
                                                          ] | tuple[
                                        tuple[sa.sage.matrix, sa.sage.matrix],
                                        tuple[sa.sage.matrix, sa.sage.matrix],
                                        tuple[sa.sage.matrix, sa.sage.matrix]]:
    """
    Reduces a regular pencil (A+tB)x = 0 to a canonical form
    [N^{u_{1}}                                     ]
    [           N^{u_{2}}                          ]
    [                     ..                       ]
    [                       ..                     ]
    [                         ..                   ]
    [                             N_{u_{s}}        ]
    [                                        J + tI]
    with N^{u} a uxu square matrix such that
    N = identity(u) + t*upper_shift_matrix(u),
    J Jordan matrix,
    I identity matrix.
    and returns
    (identity(u), upper_shift_matrix(u)), (J, I).

    If transformation is toggled on, it returns the tuple
    (P, Q), (identity(u), upper_shift_matrix(u)), (J, I)
    with (P, Q) such that
    P^{-1}*A*Q = identity(u) + J,
    P^{-1}*B*Q = upper_shift_matrix(u) + I.
    """
    assert (A.nrows() == A.ncols() == B.nrows() == A.ncols()
            and not (A + sa.var('x') * B).det().is_zero())
    EMPTY_MATRIX = sa.matrix(A.base_ring(), [])
    if A.nrows() == 0:
        if not transformation:
            return ((EMPTY_MATRIX, EMPTY_MATRIX),
                    (EMPTY_MATRIX, EMPTY_MATRIX))
        else:
            return ((EMPTY_MATRIX, EMPTY_MATRIX),
                    (EMPTY_MATRIX, EMPTY_MATRIX),
                    (EMPTY_MATRIX, EMPTY_MATRIX))
    while True:
        c = randint(-10, 10)
        if (A + c*B).det().is_zero():
            continue
        A_1 = A + c*B
        J, P_1 = (A_1.inverse() * B).jordan_form(transformation=transformation)

        J_1_LENGTH = 0
        for i in range(J.nrows() - 1, -1, -1):
            if J[i, i] != 0:
                J_1_LENGTH = i + 1
                break

        PERM = sa.block_matrix(J.base_ring(), [
            [0, sa.identity_matrix(J_1_LENGTH)],
            [sa.identity_matrix(
                J.nrows() - J_1_LENGTH), 0]])

        J_1 = J.submatrix(0, 0, J_1_LENGTH, J_1_LENGTH)
        J_0 = J.submatrix(J_1_LENGTH, J_1_LENGTH)

        if (J_1.det()).is_zero():
            continue
        break
    if J_0.ncols() <= 0 or J_0.nrows() <= 0:
        J_0 = EMPTY_MATRIX

    # J_0 may be empty, the jordan form of an empty matrix is not defined
    if J_0.nrows() > 0 and J_0.ncols() > 0:
        H, P_2 = ((sa.identity_matrix(J_0.nrows())
                   - c * J_0).inverse() * J_0).jordan_form(
                        transformation=transformation)
    else:
        H, P_2 = EMPTY_MATRIX, EMPTY_MATRIX
    J_0_IDENTITIES = sa.identity_matrix(J_0.nrows())

    if J_1.nrows() > 0 and J_1.ncols() > 0:
        M, P_3 = (J_1.inverse()
                  - c * sa.identity_matrix(J_1.nrows())).jordan_form(
                        transformation=transformation)
    else:
        M, P_3 = EMPTY_MATRIX, EMPTY_MATRIX
    J_1_IDENTITIES = sa.identity_matrix(M.nrows())

    if not transformation:
        return ((J_0_IDENTITIES, H), (M, J_1_IDENTITIES))
    else:
        L = (A_1 * P_1 *
             PERM.transpose().inverse() *
             sa.block_diagonal_matrix([P_2, P_3]))
        R = (P_1 * PERM *
             sa.block_diagonal_matrix(
                    [(sa.identity_matrix(J_0.nrows()) - c * J_0).inverse()
                     * P_2, J_1.inverse()
                     * P_3]))
        return ((L, R), (J_0_IDENTITIES, H), (M, J_1_IDENTITIES))


def kronecker_canonical_form(A: sa.sage.matrix,
                             B: sa.sage.matrix,
                             transformation=False) -> tuple[
                                 sa.sage.matrix, sa.sage.matrix] | tuple[
                                     tuple[sa.sage.matrix, sa.sage.matrix],
                                     tuple[sa.sage.matrix]]:
    """
    Computes Kronecker's canonical form of the given pencil
    (A + tB)x = 0 and return the tuple
    (KCF_A, KCF_B) such that
    KCF_A is Kronecker's canonical form of A,
    KCF_B (is Kronecker's canonical form) of B.

    If transformation is toggled on,
    it returns the tuple (P, Q), (KCF_A, KCF_B) such that
    P*A*Q = KCF_A,
    P*B*Q = KCF_B.
    """
    assert A.nrows() == B.nrows() and A.ncols() == B.ncols()

    EMPTY_MATRIX = sa.matrix(A.base_ring(), [])
    kronecker_blocks = []
    dependent_rows = []
    dependent_columns = []
    to_be_transposed = False
    L = sa.identity_matrix(A.nrows())
    R = sa.identity_matrix(A.ncols())

    if A.is_zero() and B.is_zero():
        if transformation:
            return (L, R), (A, B)
        else:
            return A, B

    A_tilde, B_tilde = A, B
    while True:
        if (A_tilde.nrows() == A_tilde.ncols() and
                not (A_tilde + sa.var('x') * B_tilde).det().is_zero()):
            break
        elif A_tilde.ncols() < A_tilde.nrows():
            to_be_transposed = True
            A_tilde = A_tilde.transpose()
            B_tilde = B_tilde.transpose()
            L, R = R.H, L.H
        (P, Q), (L_A, A_STAR), (L_B, B_STAR) = _reduction_theorem(A_tilde, B_tilde, True)
        # Check if a linear relation with constant coefficients
        # amongst the columns of the pencils has been found
        if A_STAR.nrows() == A_tilde.nrows():
            assert (L_A.is_zero() and L_B.is_zero() and
                    (P * A_tilde * Q
                     - sa.block_matrix([[L_A, A_STAR]])).is_zero() and
                    (P * B_tilde * Q
                     - sa.block_matrix([[L_B, B_STAR]])).is_zero())
            if to_be_transposed:
                dependent_rows.append((L_A.transpose(), L_B.transpose()))
            else:
                dependent_columns.append((L_A, L_B))
        else:
            assert (P * A_tilde * Q -
                    sa.block_diagonal_matrix([L_A, A_STAR])).is_zero() and (
                        P * B_tilde * Q
                        - sa.block_diagonal_matrix([L_B, B_STAR])).is_zero()
        L = sa.block_diagonal_matrix(
            [sa.identity_matrix(
                L.nrows() - P.nrows()), P]) * L
        R = R * sa.block_diagonal_matrix(
            [sa.identity_matrix(
                R.nrows() - Q.nrows()), Q])
        if to_be_transposed:
            A, B = A.transpose(), B.transpose()
            L_A = L_A.transpose()
            L_B = L_B.transpose()
            A_STAR, B_STAR = A_STAR.transpose(), B_STAR.transpose()
            L, R = R.H, L.H
            to_be_transposed = False
        if not (L_A.is_zero() and L_B.is_zero()):
            kronecker_blocks.append((L_A, L_B))
        A_tilde = A_STAR
        B_tilde = B_STAR

    P, Q = EMPTY_MATRIX, EMPTY_MATRIX
    if A_tilde.ncols() != 0 and A_tilde.nrows() != 0:
        ((P, Q), (E_1, H), (J, E_2)) = _reduce_regular_pencil(A, B, True)
        kronecker_blocks.append((E_1, H))
        kronecker_blocks.append((J, E_2))

    KCF_A = sa.block_diagonal_matrix([block[0] for block in kronecker_blocks])
    KCF_B = sa.block_diagonal_matrix([block[1] for block in kronecker_blocks])

    for row_A, row_B in dependent_rows:
        KCF_A = row_A.stack(KCF_A)
        KCF_B = row_B.stack(KCF_B)

    for col_A, col_B in dependent_columns:
        KCF_A = col_A.augment(KCF_A)
        KCF_B = col_B.augment(KCF_B)

    L = sa.block_diagonal_matrix(
        sa.identity_matrix(L.nrows() - P.nrows()), P.inverse()) * L
    R = R * sa.block_diagonal_matrix(
        sa.identity_matrix(R.nrows() - Q.nrows()), Q)

    if not transformation:
        return (KCF_A, KCF_B)
    else:
        return ((L, R), (KCF_A, KCF_B))


def stringify_pencil(A: sa.sage.matrix,
                     B: sa.sage.matrix,
                     space=30) -> str:
    """
    Returns a string representation of the given pencil.
    """
    assert A.nrows() == B.nrows() and A.ncols() == B.ncols()
    s = ""
    str_A, str_B = str(A).splitlines(), str(B).splitlines()
    s += "A:" + " " * (space+len(str_A[0])-2) + "B:\n"
    for i in range(len(str_A)):
        s += str_A[i] + " " * space + str_B[i] + "\n"
    return s


def debug() -> None:
    # Starting point is a pencil of the form (A + tB)x = 0.
    A = sa.matrix(sa.SR, [[0, 1, 0, 0], [0, 0, 1, 0]]).transpose()
    B = sa.matrix(sa.SR, [[1, 0, 0, 0], [0, 1, 0, 0]]).transpose()
    while True:
        D = sa.random_matrix(sa.ZZ, A.nrows(), A.nrows()).change_ring(sa.SR)
        if not (D.det().is_zero()):
            while True:
                C = sa.random_matrix(sa.ZZ, A.ncols(), A.ncols()).change_ring(sa.SR)
                if not (C.is_zero()):
                    break
            break
    A = D.inverse() * A * C
    B = D.inverse() * B * C
    print(f'Matrix space parent of A: {A.parent()}\n{A}\n')
    print(f'Matrix space parent of B: {B.parent()}\n{B}\n')
    (L, R), (KCF_A, KCF_B) = kronecker_canonical_form(A, B, True)

    print(stringify_pencil(KCF_A, KCF_B))
    print(stringify_pencil(L*A*R, L*B*R))
    assert (L * A * R - KCF_A).is_zero()
    assert (L * B * R - KCF_B).is_zero()


if __name__ == "__main__":
    debug()
