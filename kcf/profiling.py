import kcf_sage as kcf
import sage.all as sa
from random import randbytes, randint
import cProfile, pstats


def random_pencil(m: int, n: int) -> tuple:
    return sa.random_matrix(sa.ZZ, m, n).change_ring(sa.SR), sa.random_matrix(sa.ZZ, m, n).change_ring(sa.SR)


def batch(filename:str, size: int) -> list[tuple]:
    pr = cProfile.Profile()
    for i in range(size):
        m = randint(1, 20)
        n = randint(1, 20)
        name = f'{filename}-{i}.log'
        pr.enable()
        kcf.kronecker_canonical_form(*random_pencil(m, n))
        pr.disable()
        pstats.Stats(pr).dump_stats(name)
        with open(name, "a") as f:
            f.write(f"\nDimensions: ({m}, {n})")


def main() -> None:
    for i in range(10):
        batch(f'./profiles/{randbytes(10)}-stats-{i}', 10)


if __name__ == "__main__":
    main()