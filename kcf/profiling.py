import kcf_sage as kcf
import sage.all as sa
from random import randint
import cProfile
import pstats


def random_pencil(m: int, n: int) -> tuple:
    return (sa.random_matrix(sa.ZZ, m, n).change_ring(sa.SR),
            sa.random_matrix(sa.ZZ, m, n).change_ring(sa.SR))


def run() -> tuple[int, int, pstats.Stats]:
    pr = cProfile.Profile()
    m = randint(1, 20)
    n = randint(1, 20)
    pr.enable()
    kcf.kronecker_canonical_form(*random_pencil(m, n))
    pr.disable()
    return m, n, pstats.Stats(pr)


def main() -> None:
    for _ in range(10):
        m, n, s = run()
        print(f'Dimensions: ({m}, {n})')
        s.print_stats('kcf')


if __name__ == "__main__":
    main()
