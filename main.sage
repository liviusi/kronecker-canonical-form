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
    print(resulting_kernel)
    result = [[] for _ in range(resulting_kernel.ncols() // (degree + 1))]
    for i in range(resulting_kernel.nrows()):
        for j in range(resulting_kernel.ncols()):
            result[j // (resulting_kernel.ncols() // (degree + 1))].append(resulting_kernel[i, j])
    return (degree, sa.matrix(A.base_ring(), result).transpose())


def reduction_theorem(
        A: sa.sage.matrix,
        B: sa.sage.matrix):
    """
    Reduces a pencil (A + tB)x = 0 to a canonical form:
    [L * ]
    [0 * ].
    """
    degree, vector = compute_lowest_degree_polynomial(A, B)
    if degree <= 0:
        print("The degree of the polynomial of minimum degree " +
              "in the pencil must be greater than zero.")
        exit(1)
    print("Lowest degree polynomial:")
    print(vector)
    print("------------------------------------------------")
    print("Q:")
    Q = complete_to_a_basis(vector)
    print(Q)


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
    print(f'Matrix space parent of A: {A.parent()}')
    print(f'Matrix space parent of B: {B.parent()}')

    reduction_theorem(A, B)


if __name__ == "__main__":
    main()