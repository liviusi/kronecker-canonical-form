import sys
import sage.all as sa
from ast import literal_eval
from random import randint


RING = sa.SR
DEFAULT_PENCIL_FILE = "./pencil.txt"
EMPTY_MATRIX = sa.matrix(RING, [])


def recover_matrix(
        ring: sa.sage.rings,
        M: str) -> sa.sage.matrix:
    """
    Given a ring, constructs a matrix from its string representation M.
    """
    return sa.matrix(ring, literal_eval(M))


def recover_matrices(
        ring: sa.sage.rings,
        filename: str) -> tuple[sa.sage.matrix, sa.sage.matrix]:
    """
    Recovers the pair (A, B) from a given file.
    """
    with open(filename, "r") as f:
        lines = f.readlines()
        assert len(lines) == 2
        return (recover_matrix(ring, lines[0].split("A = ")[1]),
                recover_matrix(ring, lines[1].split("B = ")[1]))


def complete_to_a_basis(M: sa.sage.matrix) -> sa.sage.matrix:
    """
    Completes the given matrix to a basis by augmenting
    it with the identity matrix and recovering the pivots.
    """
    identity_mat = sa.matrix(M.base_ring(), [
        [1 if i == j else 0 for i in range(M.nrows())]
        for j in range(M.nrows())])
    new_M = M.augment(identity_mat)
    return new_M[:, new_M.pivots()]


def build_coefficient_matrix(
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
    of size (epsilon + 2 x epsilon + 1).
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


def compute_lowest_degree_polynomial(
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
        coefficient_matrix = build_coefficient_matrix(A, B, degree + 1)
        coefficient_matrix_kernel = (
            coefficient_matrix.right_kernel().basis_matrix()
        )
        if coefficient_matrix_kernel:
            coefficient_matrix_kernel = coefficient_matrix_kernel
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


def reduction_theorem(
        A: sa.sage.matrix,
        B: sa.sage.matrix) -> tuple[tuple[
                                            sa.sage.matrix,
                                            sa.sage.matrix, sa.sage.matrix],
                                    tuple[
                                            sa.sage.matrix,
                                            sa.sage.matrix, sa.sage.matrix]]:
    """
    Reduces a pencil (A + tB)x = 0 to a canonical form:
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
    """
    assert A.nrows() == B.nrows() and A.ncols() == B.ncols()
    degree, polynomial = compute_lowest_degree_polynomial(A, -B)
    if degree <= 0:
        print("The degree of the polynomial of minimum degree " +
              "in the pencil must be greater than zero.")
        exit(1)

    V = sa.matrix(polynomial[0].base_ring(), polynomial).transpose()
    # print(f'V:\n{V}\n')
    Q = complete_to_a_basis(V)
    P = complete_to_a_basis(A * V)
    A_tilde = P**-1 * A * Q
    B_tilde = P**-1 * B * Q

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
        return ((L_A, EMPTY_MATRIX, EMPTY_MATRIX),
                (L_B, EMPTY_MATRIX, EMPTY_MATRIX))

    rows = []
    for i in range(L_A.nrows()):
        row = []
        for j in range(L_A.ncols(), A.ncols()):
            row.append(A_tilde[i, j])
        rows.append(row)
    D = sa.matrix(A.base_ring(), rows)

    rows = []
    for i in range(L_A.nrows(), A.nrows()):
        row = []
        for j in range(L_A.ncols(), A.ncols()):
            row.append(A_tilde[i, j])
        rows.append(row)
    A_STAR = sa.matrix(A.base_ring(), rows)

    rows = []
    for i in range(L_B.nrows()):
        row = []
        for j in range(L_B.ncols(), B.ncols()):
            row.append(B_tilde[i, j])
        rows.append(row)
    F = sa.matrix(B.base_ring(), rows)

    rows = []
    for i in range(L_B.nrows(), B.nrows()):
        row = []
        for j in range(L_B.ncols(), B.ncols()):
            row.append(B_tilde[i, j])
        rows.append(row)
    B_STAR = sa.matrix(B.base_ring(), rows)

    W = Z = Y = None

    if degree - 1 > 0:
        M = build_coefficient_matrix(A_STAR, B_STAR, degree-1)
        if degree > 1:
            rows = []
            sign = -1
            for i in range(F.nrows() - 1):
                row = []
                sign *= -1
                for j in range(F.ncols()):
                    row.append(sign * (F[i+1, j] - D[i, j]))
                rows.append(row)
            W = sa.matrix(A.base_ring(), rows)
            Z = M.solve_left(W)
            assert Z.ncols() % degree == 0
            rows = []
            sign = -1
            rounds = 0
            for i in range(degree):
                row = []
                sign *= -1
                for j in range(Z.ncols() // degree):
                    row.append(sign * Z[0, min(i + j + len(rows), Z.ncols()-1)])
                rows.append(row)
                rounds += 1
            Y = sa.matrix(A.base_ring(), rows)
    else:
        Y = sa.matrix(A.base_ring(), [[1 for _ in range(B_STAR.nrows())]])

    X = [[None for _ in range(A_STAR.ncols())] for _ in range(degree + 1)]
    for i in range(degree):
        for j in range(A_STAR.ncols()):
            X[i][j] = F[i, j] + Y.row(i) * B_STAR.column(j)
    for j in range(A_STAR.ncols()):
        X[degree][j] = D[degree - 1, j] + Y.row(degree - 1) * A_STAR.column(j)

    X = sa.matrix(A.base_ring(), X)
    return ((L_A, D, A_STAR), (L_B, F, B_STAR))
"""
    print(f'L_A:\n{L_A}\nA_STAR:\n{A_STAR}\nB_STAR:\n{B_STAR}\nD:\n{D}\nF:\n{F}\nY:\n{Y}\nX:\n{X}\n')
    print(f'Y*A_STAR:\n{Y*A_STAR}\nD-L_A*X:\n{D-L_A*X}')
    print(f'RESULT:\n{D + Y * A_STAR - L_A * X}\n')

    assert ((D + Y * A_STAR - L_A * X).is_zero()) and (
        (F + Y * B_STAR - L_B * X).is_zero())
"""


def reduce_regular_pencil(A: sa.sage.matrix,
                          B: sa.sage.matrix) -> tuple[
                                        tuple[sa.sage.matrix, sa.sage.matrix],
                                        tuple[sa.sage.matrix, sa.sage.matrix]]:
    assert A.nrows() == A.ncols() == B.nrows() == A.ncols()
    if A.nrows() == 0:
        return ((EMPTY_MATRIX, EMPTY_MATRIX), (EMPTY_MATRIX, EMPTY_MATRIX))
    while True:
        c = randint(-10, 10)
        if (A + c*B).det().is_zero():
            continue
        A_1 = A + c*B
        J = (A_1.inverse() * B).jordan_form()

        J_1_LENGTH = -1
        for i in range(J.nrows() - 1, -1, -1):
            if J[i, i] != 0:
                J_1_LENGTH = i + 1
                break

        J_1 = J.submatrix(0, 0, J_1_LENGTH, J_1_LENGTH)
        J_0 = J.submatrix(J_1_LENGTH, J_1_LENGTH,
                          J.nrows()-J_1.nrows(), J.ncols()-J_1.ncols())

        H = ((sa.identity_matrix(J_0.nrows())
              - c * J_0).inverse() * J_0).jordan_form() if J_0.nrows() > 0 else EMPTY_MATRIX
        J_0_IDENTITIES = sa.identity_matrix(J_0.nrows())

        M = (J_1.inverse() - c * sa.identity_matrix(J_1.nrows())).jordan_form()
        J_1_IDENTITIES = sa.identity_matrix(M.nrows())

        return ((J_0_IDENTITIES, H), (M, J_1_IDENTITIES))


def kcf(A: sa.sage.matrix, B: sa.sage.matrix):
    """
    Computes the kronecker canonical form of the given pencil
    (A + tB)x = 0.
    """
    kronecker_blocks = []
    to_be_transposed = False
    KCF = None
    while True:
        if (A.nrows() == A.ncols() and
                not (A + sa.var('x') * B).det().is_zero()):
            break
        elif A.ncols() < A.nrows():
            to_be_transposed = True
            A = A.transpose()
            B = B.transpose()
        (L_A, _, A_STAR), (L_B, _, B_STAR) = reduction_theorem(A, B)
        L = L_A + sa.var('t') * L_B
        if to_be_transposed:
            L = L.transpose()
            to_be_transposed = False
        kronecker_blocks.append(L)
        A = A_STAR
        B = B_STAR

    if A.ncols() != 0 and A.nrows() != 0:
        ((E_1, H), (J, E_2)) = reduce_regular_pencil(A, B)
        kronecker_blocks.append(E_1 + sa.var('t') * H)
        kronecker_blocks.append(J + sa.var('t') * E_2)

    KCF = sa.block_diagonal_matrix(kronecker_blocks)
    return KCF


def main() -> None:
    """
    Main function.
    """
    filename = None
    if (len(sys.argv) <= 1):
        print("No input file has been selected. Default will be used:" +
              f" {DEFAULT_PENCIL_FILE}.")
        filename = DEFAULT_PENCIL_FILE
    elif len(sys.argv) != 2:
        print("Only a single argument may be provided." +
              " Usage: sage main.sage <pencil_filename>")
        exit(1)
    else:
        filename = sys.argv[1]

    # Starting point is a pencil of the form (A + tB)x = 0.
    # TODO: this crashes because of incompatible matrices sizes in reduction_theorem()
    # A, B = sa.random_matrix(sa.ZZ, 10, 8).change_ring(RING), sa.random_matrix(sa.ZZ, 10, 8).change_ring(RING)
    A, B = sa.random_matrix(sa.ZZ, 9, 6).change_ring(RING), sa.random_matrix(sa.ZZ, 9, 6).change_ring(RING)

    # Matrices sizes must be compatible.
    assert A.nrows() == B.nrows() and B.ncols() == A.ncols()

    # Info:
    print(f'Matrix space parent of A: {A.parent()}\n{A}\n')
    print(f'Matrix space parent of B: {B.parent()}\n{B}\n')
    KCF = kcf(A, B)
    print(f'KCF:\n{KCF}\n')


if __name__ == "__main__":
    main()
