import sys
import sage.all as sa
from ast import literal_eval
from random import randint
from typing import List


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
        B: sa.sage.matrix,
        transformation=False):
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
        M = build_coefficient_matrix(A_STAR, B_STAR, degree-1)
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
        left_M = sa.block_matrix(sa.SR, [[sa.identity_matrix(degree), Y], [0, sa.identity_matrix(A_tilde.nrows() - degree)]])
        right_M = sa.block_matrix(sa.SR, [[sa.identity_matrix(degree + 1), -X], [0, sa.identity_matrix(A_tilde.ncols() - degree - 1)]])
        return ((left_M* P**-1, Q * right_M), (L_A, A_STAR), (L_B, B_STAR))


def reduce_regular_pencil(A: sa.sage.matrix,
                          B: sa.sage.matrix,
                          transformation=False) -> tuple[
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

        J_1_LENGTH = 0
        for i in range(J.nrows() - 1, -1, -1):
            if J[i, i] != 0:
                J_1_LENGTH = i + 1
                break

        J_1 = J.submatrix(0, 0, J_1_LENGTH, J_1_LENGTH)
        J_0 = J.submatrix(J_1_LENGTH-1, J_1_LENGTH-1,
                          J.nrows()-J_1.nrows(), J.ncols()-J_1.ncols())

        # J_0 may be empty, the jordan form of an empty matrix is not defined
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
    assert A.nrows() == B.nrows() and A.ncols() == B.ncols()
    kronecker_blocks = []
    to_be_transposed = False
    L = sa.identity_matrix(A.nrows())
    R = sa.identity_matrix(A.ncols())
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
        (P, Q), (L_A, A_STAR), (L_B, B_STAR) = reduction_theorem(A_tilde, B_tilde, True)
        assert (P * A_tilde * Q - sa.block_diagonal_matrix([L_A, A_STAR])).is_zero() and (P * B_tilde *  Q - sa.block_diagonal_matrix([L_B, B_STAR])).is_zero()
        L = sa.block_diagonal_matrix([sa.identity_matrix(L.nrows() - P.nrows()), P]) * L
        R = R * sa.block_diagonal_matrix([sa.identity_matrix(R.nrows() - Q.nrows()), Q])
        if to_be_transposed:
            A, B = A.transpose(), B.transpose()
            L_A = L_A.transpose()
            L_B = L_B.transpose()
            A_STAR, B_STAR = A_STAR.transpose(), B_STAR.transpose()
            L, R = R.H, L.H
            to_be_transposed = False
        kronecker_blocks.append((L_A, L_B))
        A_tilde = A_STAR
        B_tilde = B_STAR

    if A_tilde.ncols() != 0 and A_tilde.nrows() != 0:
        ((E_1, H), (J, E_2)) = reduce_regular_pencil(A, B)
        kronecker_blocks.append((E_1, H))
        kronecker_blocks.append((J, E_2))

    KCF_A = sa.block_diagonal_matrix([block[0] for block in kronecker_blocks])
    KCF_B = sa.block_diagonal_matrix([block[1] for block in kronecker_blocks])

    # assert (L*A*R - KCF_A).is_zero() and (L*B*R - KCF_B).is_zero()

    return (KCF_A, KCF_B)

def print_pencil(A: sa.sage.matrix, B: sa.sage.matrix, space=30):
    assert A.nrows() == B.nrows() and A.ncols() == B.ncols()
    str_A, str_B = str(A).splitlines(), str(B).splitlines()
    print("A:" + " " * (space+len(str_A[0])-2) + "B:")
    for i in range(len(str_A)):
        print(str_A[i] + " " * space + str_B[i])


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
    A, B = recover_matrices(RING, filename)

    # Info:
    print(f'Matrix space parent of A: {A.parent()}\n{A}\n')
    print(f'Matrix space parent of B: {B.parent()}\n{B}\n')
    KCF_A, KCF_B = kcf(A, B)
    print_pencil(KCF_A, KCF_B)


if __name__ == "__main__":
    main()
