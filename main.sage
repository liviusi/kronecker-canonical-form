import sage.all as sa


def compute_kernel(
        A: sa.sage.matrix,
        B: sa.sage.matrix) -> tuple[int, sa.sage.matrix]:
    """
    Returns the pair (degree, v) with v representing the
    polynomial of minimum degree in the kernel of the pencil
    (A + tB)x = 0
    and degree its degree.
    """
    assert A.nrows() == B.nrows() and B.ncols() == A.ncols()
    assert A.base_ring() == B.base_ring()
    degree = 0
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
                return (degree, right_kernel_A)
        else:
            rows = []
            for i in range(degree + 1):
                rows.append([A if i == j else
                             (B if i == j + 1 else 0)
                             for j in range(degree + 1)])
            rows.append([B if j == degree else 0 for j in range(degree+1)])
            coefficient_matrix = sa.block_matrix(A.base_ring(), rows)
            print(f"{degree=}")
            # print(coefficient_matrix)
            matrix_kernel = coefficient_matrix.right_kernel().basis_matrix()
            if matrix_kernel:
                return (degree, matrix_kernel)
        degree += 1


def reduction_theorem(
        A: sa.sage.matrix,
        B: sa.sage.matrix):
    """
    Reduces a pencil (A + tB)x = 0 to a canonical form:
    [L * ]
    [0 * ].
    """
    degree, vector = compute_kernel(A, B)
    k = len(vector.list())
    if k <= 0:
        print("The degree of the polynomial of minimum degree " +
              "in the pencil must be greater than zero.")
        exit(1)
    assert not len(A.base_ring().linear_dependence(vector)) <= 0
    if k < A.ncols():
        vector = sa.vector(vector.base_ring(), vector.list() +
                           [0 for _ in range(A.ncols() - len(vector.list()))])
    w_matrix = [A.solve_right(vector[i]) for i in range(k)]
    if k < A.ncols():
        w_matrix += [sa.vector(A.sa.base_ring(),
                               [0 for _ in range(A.ncols() - k)])]
    return (vector, sa.matrix(A.base_ring(), w_matrix))


# Starting point is a pencil of the form (A + tB)x = 0.
ring = sa.AA

A = sa.matrix(ring, [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]])
B = sa.matrix(ring, [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])

assert A.nrows() == B.nrows() and B.ncols() == A.ncols()
print(compute_kernel(A, B))
