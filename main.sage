from re import I
import sage.all as sa
import sys
from ast import literal_eval


ring = sa.AA
DEFAULT_PENCIL_FILE = "./pencil.txt"


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
    Completes the given matrix to a basis.
    TODO: actually document how this function works.
    """
    identity_mat = sa.matrix(M.base_ring(), [
        [1 if i == j else 0 for i in range(M.nrows())]
        for j in range(M.nrows())])
    new_M = M.augment(identity_mat)
    # assert pivots are in the very first positions
    return new_M[:, new_M.pivots()]


def build_coefficient_matrix(
        A: sa.sage.matrix,
        B: sa.sage.matrix,
        epsilon: int, change_sign=False) -> sa.sage.matrix:
    """
    Given two matrices (A, B) of the same size returns a matrix of
    the following form:
    [A   0     0]
    [B   A      ]
    [0   B      ]
    [          A]
    [0 0  ...  B]
    of size (epsilon + 2 x epsilon + 1).
    change_sign is to be toggled on if the matrix is to be built
    using -B instead of B.
    """
    assert A.ncols() == B.ncols() and A.nrows() == B.nrows()
    if epsilon <= 0:
        return sa.matrix(A.base_ring(), [])
    else:
        sign = -1 if change_sign else 1
        rows = []
        if epsilon == 1:
            rows.append([A])
            rows.append([sign * B])
        else:
            for i in range(epsilon):
                rows.append([A if i == j else
                             (sign * B if i == j + 1 else 0)
                             for j in range(epsilon)])
            rows.append([sign * B if j == epsilon - 1 else 0 for j in range(epsilon)])
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
        coefficient_matrix = build_coefficient_matrix(A, B, degree + 1, True)
        coefficient_matrix_kernel = (
            coefficient_matrix.right_kernel().basis_matrix()
        )
        if coefficient_matrix_kernel:
            coefficient_matrix_kernel = coefficient_matrix_kernel
            break
        degree += 1
    coefficient_matrix_kernel = coefficient_matrix_kernel.transpose()
    indices = [i for i in range(0, coefficient_matrix_kernel.nrows(), coefficient_matrix_kernel.nrows() // (degree + 1))]
    coefficient_matrix_kernel.subdivide(indices, None)
    result = []
    for index in indices:
        col = []
        for i in range(index, index + coefficient_matrix_kernel.nrows() // (degree + 1)):
            col.append(coefficient_matrix_kernel[i, 0])
        result.append(sa.vector(ring, col))

    return (degree, result)


def reduction_theorem(
        A: sa.sage.matrix,
        B: sa.sage.matrix) -> tuple[tuple[
                                            sa.sage.matrix, sa.sage.matrix,
                                            sa.sage.matrix, sa.sage.matrix],
                                    tuple[
                                            sa.sage.matrix, sa.sage.matrix,
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
    degree, polynomial = compute_lowest_degree_polynomial(A, B)
    if degree <= 0:
        print("The degree of the polynomial of minimum degree " +
              "in the pencil must be greater than zero.")
        exit(1)

    V = sa.matrix(polynomial[0].base_ring(), polynomial).transpose()
    print(f'V:\n{V}\n')
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
        return ((L_A, sa.matrix(A.base_ring(), []), sa.matrix(A.base_ring(), [])),
                        (L_B, sa.matrix(B.base_ring(), []),
                                sa.matrix(B.base_ring(), []), sa.matrix(B.base_ring(), [])))

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

    # Does this really need to be computed each and every time?
    # It is already shown how such matrices can always be defined,
    # is this even useful?
    W = Z = Y = None

    if degree - 1 > 0:
        M = build_coefficient_matrix(A_STAR, B_STAR, degree-1)
        # print(f'M:\n{M}\n')
        if degree > 1:
            rows = []
            sign = -1
            for i in range(F.nrows() - 1):
                row = []
                for j in range(F.ncols()):
                    sign *= -1
                    row.append(sign * (F[i+1, j] - D[i, j]))
                rows.append(row)
            W = sa.matrix(A.base_ring(), rows)
            Z = M.solve_left(W)
            assert Z.ncols() % degree == 0
            rows = []
            for i in range(Z.ncols()):
                row = []
                sign = -1
                for k in range(Z.ncols() // degree):
                    sign *= -1
                    row.append(sign * Z[i])
                rows.append(row)
            Y = sa.matrix(A.base_ring(), rows)
    else:
        Y = sa.matrix(A.base_ring(), [[1 for _ in range(degree)]])

    X = [[None for _ in range(A_STAR.ncols())] for _ in range(degree + 1)]
    for i in range(degree):
        for j in range(A_STAR.ncols()):
            X[i][j] = F[i, j] + Y.row(i) * B_STAR.column(j)
    for j in range(A_STAR.ncols()):
        X[degree][j] = D[degree - 1, j] + Y.row(degree - 1) * A_STAR.column(j)

    X = sa.matrix(A.base_ring(), X)

    assert ((D + Y * A_STAR - L_A * X).is_zero()) and ((F + Y * B_STAR - L_B * X).is_zero())

    return ((L_A, D, A_STAR), (L_B, F, B_STAR))


"""
    rows = []
    for i in range(degree):
        row = []
        for k in range(F.ncols()):
            print(f'{F[i,k]}')
            row.append(F[i+1, k] - D[i, k])
        rows.append(row)
    W = sa.matrix(A.base_ring(), rows)
    print(W)
"""



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
    A, B = recover_matrices(ring, filename)

    # Matrices sizes must be compatible.
    assert A.nrows() == B.nrows() and B.ncols() == A.ncols()

    # Info:
    print(f'Matrix space parent of A: {A.parent()}\n{A}\n')
    print(f'Matrix space parent of B: {B.parent()}\n{B}\n')

    (L_A, D, A_STAR), (L_B, F, B_STAR) = reduction_theorem(A, B)
    print(f"L_A:\n{L_A}\n")
    print(f"D:\n{D}\n")
    print(f"A_STAR:\n{A_STAR}\n")
    print(f"L_B:\n{L_B}\n")
    print(f"F:\n{F}\n")
    print(f"B_STAR:\n{B_STAR}\n")


if __name__ == "__main__":
    main()
