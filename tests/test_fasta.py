from __future__ import annotations

import random
import sys

from Bio import SeqIO

from nnfasta import nnfastas


def _ok(r, rec):
    assert r.id == rec.id
    assert r.seq == str(rec.seq)
    assert r.description == rec.description


def _test():

    ffs = sys.argv[1:]
    if not ffs:
        print("no fasta files specified")
        return

    fasta = nnfastas(ffs)
    print("checking list of filenames")
    _check(ffs, fasta)

    # check raw bytes
    blist = [open(f, "rb").read() for f in ffs]
    fasta = nnfastas(blist)
    print("checking list of bytes")
    _check(ffs, fasta)

    # check open files
    iolist = [open(f, "rb") for f in ffs]
    fasta = nnfastas(iolist)
    print("checking list of open files")
    _check(ffs, fasta)

    print("checking to_rec")
    _check2(ffs, fasta)
    print("done.")


def _check(ffs, fasta):
    full = []
    for ff in ffs:
        with open(ff, encoding="utf8") as h:
            full.extend(SeqIO.parse(h, "fasta"))

    total = len(full)
    assert total == len(fasta)

    for i, rec in enumerate(full):
        r = fasta[i]
        _ok(r, rec)

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


def _check2(ffs: list[str], fasta):
    full = []
    for ff in ffs:
        with open(ff, encoding="utf8") as h:
            full.extend(SeqIO.parse(h, "fasta"))

    for f1, f2 in zip(full, fasta):
        s2 = f2.to_rec()
        assert f1.id == s2.id
        assert f1.seq == s2.seq
        assert f1.description == s2.description


_test()
