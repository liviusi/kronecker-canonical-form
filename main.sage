import sage.all as sa
import sys
from ast import literal_eval
from itertools import chain

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
    return new_M[:, new_M.pivots()]


def compute_lowest_degree_polynomial(
        A: sa.sage.matrix,
        B: sa.sage.matrix) -> tuple[int, list[sa.sage.matrix]]:
    """
    Returns the pair (degree, v) with v representing the
    polynomial of minimum degree in the kernel of the pencil
    (A + tB)x = 0
    and degree its degree.
    """
    assert A.nrows() == B.nrows() and B.ncols() == A.ncols()
    assert A.base_ring() == B.base_ring()
    degree = 0
    resulting_kernel = None
    while True:
        if degree == 0:
            """
            [A] v = 0
            [B]
            """
            # Computing: A * v_0 = 0
            right_kernel_A = A.right_kernel().basis_matrix()
            if not right_kernel_A:
                print("Ker(A) is empty.")
                exit(1)
            right_kernel_B = B.right_kernel().basis_matrix()
            if not right_kernel_B:
                print("Ker(B) is empty.")
                exit(1)
            if (right_kernel_A == right_kernel_B):
                resulting_kernel = right_kernel_A
                break
        else:
            rows = []
            for i in range(degree + 1):
                rows.append([A if i == j else
                             (B if i == j + 1 else 0)
                             for j in range(degree + 1)])
            rows.append([B if j == degree else 0 for j in range(degree+1)])
            coefficient_matrix = sa.block_matrix(A.base_ring(), rows)
            coefficient_matrix_kernel = (
                coefficient_matrix.right_kernel().basis_matrix()
            )
            if coefficient_matrix_kernel:
                resulting_kernel = coefficient_matrix_kernel
                break
        degree += 1
    resulting_kernel = resulting_kernel.transpose()
    indices = [i for i in range(0, resulting_kernel.nrows(), resulting_kernel.nrows() // (degree + 1))]
    resulting_kernel.subdivide(indices, None)
    result = []
    for index in indices:
        v = []
        for j in range(resulting_kernel.ncols()):
            col = []
            for i in range(index, index + resulting_kernel.nrows() // (degree + 1)):
                col.append(resulting_kernel[i, j])
            v.append(col)
        v = sa.matrix(ring, v).transpose()
        result.append(v)

    return (degree, result)


def reduction_theorem(
        A: sa.sage.matrix,
        B: sa.sage.matrix) -> tuple[
                                        sa.sage.matrix, sa.sage.matrix,
                                        sa.sage.matrix, sa.sage.matrix]:
    """
    Reduces a pencil (A + tB)x = 0 to a canonical form:
    A_tilde = P^-1 * A * Q =
    [0 1 0 ... 0 |        ]
    [0 0 1 ... 0 |   D    ]
    [0 0 0 ... 1 |        ]
    [------------|------- ]
    [     0      | A_star ]

    B_tilde = P^-1 * B * Q =
    [1 0 0 ... 0 |        ]
    [0 1 0 ... 0 |   F    ]
    [0 0 ... 1 0 |        ]
    [------------|------- ]
    [     0      | B_star ]
    """
    degree, polynomial = compute_lowest_degree_polynomial(A, B)
    if degree <= 0:
        print("The degree of the polynomial of minimum degree " +
              "in the pencil must be greater than zero.")
        exit(1)

    V = sa.block_matrix(ring, [[v for v in polynomial]], subdivide=True)
    print(f'V:\n{V}\n')
    print(f'Product:\n{A*polynomial[0]}\n')
    print(f'Product:\n{B*polynomial[len(polynomial)-1]}\n')
    Q = complete_to_a_basis(V)
    P = complete_to_a_basis(A * V)
    A_tilde = P ^ -1 * A * Q
    B_tilde = P ^ -1 * B * Q
    return (A_tilde, B_tilde, P, Q)


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

    A_tilde, B_tilde, P, Q = reduction_theorem(A, B)
    print(f"A_tilde:\n{A_tilde}\n")
    print(f"B_tilde:\n{B_tilde}\n")
    print(f"P:\n{P}\n")
    print(f"Q:\n{Q}\n")


if __name__ == "__main__":
    main()
