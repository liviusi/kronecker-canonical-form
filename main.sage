import sage.all as sa


def compute_kernel(
        A: sa.sage.matrix,
        B: sa.sage.matrix):
    """
    Returns the vector representing the polynomial of minimum degree
    in the kernel of the pencil
    (A + tB)x = 0.
    """
    assert A.nrows() == B.nrows() and B.ncols() == A.ncols()
    found_root = False
    vectors = list()
    while not found_root:
        if len(vectors) == 0:
            # Computing: A * v_0 = 0
            tmp = [item for item in list(A.right_kernel().basis_matrix())[0]]
            if not tmp:
                print("Ker(A) is empty.")
                exit(1)
            right_kernel_A = sa.vector(A.base_ring(), tmp)
            vectors.append(right_kernel_A)
            tmp = [item for item in list(B.right_kernel().basis_matrix())[0]]
            if not tmp:
                print("Ker(B) is empty.")
                exit(1)
            right_kernel_B = sa.vector(B.base_ring(), tmp)
            if (right_kernel_A == right_kernel_B):
                found_root = True
        else:
            # Computing: A * v_i + B * v_{i-1} = 0
            # [A * v_i | B * v_{i-1}] y = 0 has a solution if and only if
            # the polynomial is of this degree.
            break
    return vectors


# Starting point is a pencil of the form (A + tB)x = 0.
ring = sa.AA

A = sa.matrix(ring, [[0, 1, 2], [2, 3, 6], [4, 2, 4], [0, 2, 4]])
B = sa.matrix(ring, [[1, 4, 6], [2, 8, 5], [2, 8, 12], [1, 0, 1]])

assert A.nrows() == B.nrows() and B.ncols() == A.ncols()
