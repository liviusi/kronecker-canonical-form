import sage.all as sa


def remove_linear_dependency(
        M: sa.sage.matrix,
        ring: sa.sage.rings,
        remove_rows=True,
        remove_columns=False) -> sa.sage.matrix:
    """
    Returns the matrix M deprived of its linear dependent rows
    and or columns if any.
    """
    resulting_matrix = M
    while True:
        new_items = []
        removed_dependency = False
        if remove_rows:
            V = ring ^ resulting_matrix.ncols()
            for row in resulting_matrix.rows():
                r = sa.vector(ring, [item for item in list(row)])
                tmp_rows = new_items.copy()
                tmp_rows.append(r)
                if (not new_items or
                        (not V.linear_dependence(tmp_rows))):
                    new_items.append(r)
                else:
                    removed_dependency = True
            if removed_dependency:
                resulting_matrix = sa.matrix(ring, new_items)
                print(resulting_matrix)
                continue
            else:
                new_items = []
        if remove_columns:
            V = ring ^ resulting_matrix.nrows()
            for column in resulting_matrix.columns():
                c = sa.vector(ring, [item for item in list(column)])
                tmp_columns = new_items.copy()
                tmp_columns.append(c)
                if (not new_items or
                        (not V.linear_dependence(tmp_columns))):
                    new_items.append(c)
                else:
                    removed_dependency = True
            if removed_dependency:
                new_rows = [[row[i] for row in new_items] for i
                            in range(len(new_items[0]))]
                resulting_matrix = sa.matrix(ring, new_rows)
                continue
        break
    return resulting_matrix


# TODO: fix this
# ValueError: number of rows of self must equal degree of right-hand side.
def compute_degree(
        A: sa.sage.matrix,
        B: sa.sage.matrix,
        ring: sa.sage.rings) -> int:
    """
    Returns the minimum degree of the vector of polynomials in the kernel
    of the pencil (A + tB)x = 0.
    """
    assert A.nrows() == B.nrows() and B.ncols() == A.ncols()
    found_root = False
    vectors = list()
    tmp_B = remove_linear_dependency(B, ring)
    print(tmp_B)
    while not found_root:
        if len(vectors) == 0:
            # Computing: A * v_0 = 0
            tmp = [item for item in list(A.right_kernel().basis_matrix())[0]]
            if not tmp:
                print("Ker(A) is empty.")
                exit(1)
            vectors.append(sa.vector(ring, tmp))
        else:
            # Computing: A * v_i = -B * v_{i-1}
            v = A.solve_right(-(tmp_B*vectors[-1]))
            vectors.append(v)
        print(vectors[-1])
        if (vectors[-1] == tmp_B.solve_right(vectors[-1])):
            print(vectors)
            found_root = True
    return len(vectors) - 1


# Starting point is a pencil of the form (A + tB)x = 0.
ring = sa.RR
A = sa.matrix(ring, [[0, 1, 2], [2, 3, 6], [4, 2, 4], [0, 2, 4]])
B = sa.matrix(ring, [[1, 4, 6], [2, 8, 5], [2, 8, 12], [1, 0, 1]])

assert A.nrows() == B.nrows() and B.ncols() == A.ncols()

degree = compute_degree(A, B, ring)
print(f"{degree=}")
