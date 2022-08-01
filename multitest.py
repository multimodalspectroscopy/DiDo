import multiprocessing as mp
import time
import os
from unittest import result


def sum_up_to(number):
    return sum(range(1, number + 1))

tic = time.perf_counter()

array = [1,2,3]
def f(x):
    return x*x


if __name__ == '__main__':
    N= mp.cpu_count()
    print(N)
    with mp.Pool(processes = N) as p:
        results = p.map(f, range(int(10e1)))
toc = time.perf_counter()
print(toc-tic)

