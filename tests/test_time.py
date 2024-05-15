from __future__ import annotations

import sys
import timeit

from Bio import SeqIO

from nnfasta import nnfastas


def _test_bio():
    fastas = sys.argv[1:]
    n = 0
    for fasta in fastas:
        with open(fasta, encoding="utf8") as fp:
            for rec in SeqIO.parse(fp, "fasta"):
                n += len(rec.seq)

    # print(n)


def _test_nn():
    fastas = sys.argv[1:]
    n = 0
    for rec in nnfastas(fastas):
        n += len(rec.seq)


def _test():
    # seems like nnfasta is 10x slower
    print("Bio:", timeit.timeit("_test_bio()", globals=globals(), number=10))
    print("nnfasta:", timeit.timeit("_test_nn()", globals=globals(), number=10))


_test()
