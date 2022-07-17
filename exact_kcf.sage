import sage.all as sa


def cleanup_matrix(
        M: sa.sage.matrix,
        epsilon) -> sa.sage.matrix:
    """
    Finds and replaces all elements lesser than epsilon and replaces
    all of them with zero.
    """
    rows = []
    for row in M.rows():
        rows.append([0 if abs(r) <= epsilon else r for r in row])
    return sa.matrix(M.base_ring(), rows)


def compress_columns(
        M: sa.sage.matrix,
        epsilon) -> tuple[int, sa.sage.matrix, sa.sage.matrix]:
    """
    Returns rank of a matrix M and the factors X, V of a decomposition
    of the following type:
    M * V = [X | 0]
    X has rank(M) columns.
    """
    _, SIGMA, V = M.SVD()
    tmp = cleanup_matrix(M * V, epsilon)
    rank = sum(map(lambda x: abs(x) > epsilon, SIGMA.diagonal()))
    rows = []
    for row in tmp.rows():
        rows.append(row[0:rank])
    return (rank, sa.matrix(M.base_ring(), rows), cleanup_matrix(V, epsilon))


def compress_rows(
        M: sa.sage.matrix,
        epsilon) -> tuple[int, sa.sage.matrix, sa.sage.matrix]:
    """
    Returns the rank of a matrix M and the factors
    U_STAR, X of a decomposition of the following type:
    U_STAR * M = [ X ]
                 [ 0 ]
    X has rank(M) rows.
    """
    U, SIGMA, _ = M.SVD()
    U_STAR = cleanup_matrix(U.conjugate_transpose(), epsilon)
    tmp = cleanup_matrix(U_STAR * M, epsilon)
    rank = sum(map(lambda x: abs(x) > epsilon, SIGMA.diagonal()))
    rows = []
    for i in range(rank):
        rows.append(tmp[i])
    return (rank, U_STAR, sa.matrix(M.base_ring(), rows))


# TODO: fix this
def exact_kcf_column_structure(
        A: sa.sage.matrix,
        B: sa.sage.matrix,
        epsilon):
    """
    Computes the Kronecker column structure of the matrix pencil
    (A + tB)x = 0.
    """
    assert A.nrows() == B.nrows() and A.ncols() == B.ncols()
    assert A.nrows() != 0 and A.ncols() != 0
    s = []
    r = []
    new_A = A
    new_B = B
    for i in range(A.ncols()):
        m, n = A.nrows(), A.ncols()
        if n == 0 or m == 0:
            break
        rank_B, X_B, _ = compress_columns(B, epsilon)
        s.append(n - rank_B)
        if rank_B == n:
            break
        # Transform A
        new_A *= X_B
        # Transform B
        new_B *= X_B
        # Row compress the part of A in the kernel of B
        rank_A, X_A, _ = compress_rows(A, epsilon)
        r.append(rank_A)
        new_A = cleanup_matrix(A * X_A, epsilon)
        new_B = cleanup_matrix(B * X_A, epsilon)
    return (new_A, new_B)


# Starting point is a pencil of the form (A + tB)x = 0.
ring = sa.CDF
A = sa.matrix(ring, [[0, 1, 2], [2, 3, 6], [4, 2, 4], [0, 2, 4]])
B = sa.matrix(ring, [[1, 4, 6], [2, 8, 5], [2, 8, 12], [1, 0, 1]])

assert A.nrows() == B.nrows() and B.ncols() == A.ncols()
# print(exact_kcf_column_structure(A, B, pow(10, -11)))
