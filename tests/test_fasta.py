import sys
import random
from nnfasta import nnfastas
from Bio import SeqIO


def _ok(r, rec):
    assert r.id == rec.id
    assert r.seq == str(rec.seq)
    assert r.description == rec.description


def _test():

    ffs = sys.argv[1:]
    if not ffs:
        print("no fasta files")
        return

    fasta = nnfastas(ffs)

    full = []
    for ff in ffs:
        with open(ff, encoding="utf8") as h:
            full.extend(SeqIO.parse(h, "fasta"))

    assert len(full) == len(fasta)

    for i, rec in enumerate(full):
        r = fasta[i]
        _ok(r, rec)

    total = len(full)
    for s, e in [
        (s, random.randint(s + 1, total))
        for s in random.sample(range(0, total), min(2000, total))
    ]:
        if s == e:
            continue
        recs = full[s:e]
        rs = fasta[s:e]
        for rec, r in zip(recs, rs):
            _ok(r, rec)

    for rl in [
        random.sample(range(0, total), random.randint(10, 60)) for _ in range(40)
    ]:
        rs = fasta[rl]
        recs = [full[i] for i in rl]
        for rec, r in zip(recs, rs):
            _ok(r, rec)
    print(f"OK {total}")


_test()
