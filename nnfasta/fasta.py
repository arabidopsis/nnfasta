import mmap
import re
import os
import bisect
from dataclasses import dataclass
from collections.abc import Sequence
from typing import Iterator, overload


@dataclass
class Record:
    id: str
    description: str
    seq: str

    @property
    def name(self) -> str:
        return self.id


PREFIX = re.compile(b"(\n>|\r>|^>)", re.M)

WHITE = re.compile(rb"\W+", re.I | re.M)
EOL = re.compile(rb"\n|\r", re.I | re.M)


def remove_white(s: bytes) -> bytes:
    return WHITE.sub(b"", s)


def nnfastas(
    fasta_files: Sequence[os.PathLike], encoding: str | None = None
) -> Sequence[Record]:
    if not fasta_files:
        raise ValueError("no fasta files!")
    if len(fasta_files) == 1:
        return RandomFasta(fasta_files[0], encoding=encoding)
    return CollectionFasta(fasta_files, encoding=encoding)


class RandomFasta(Sequence[Record]):
    """Memory mapped fasta file."""

    ENCODING = "ascii"

    def __init__(self, fasta_file: os.PathLike, encoding: str | None = None):

        self.encoding = encoding or self.ENCODING
        self.fp = open(fasta_file, "rb")  # pylint: disable=consider-using-with
        self.fasta = mmap.mmap(self.fp.fileno(), 0, prot=mmap.PROT_READ)
        self.pos = self._find_pos()

    def __del__(self):
        if self is not None and self.fp:
            self.fp.close()
            self.fp = None

    def _find_pos(self) -> list[tuple[int, int]]:
        f = [(h.start(), h.end()) for h in PREFIX.finditer(self.fasta)]
        end, start = zip(*f)
        end = end[1:] + (len(self.fasta),)
        return list(zip(start, end))

    def get_idx(self, idx: int) -> Record:
        s, e = self.pos[idx]
        b = self.fasta[s:e]  # mmap go to disk
        m = EOL.search(b)
        if not m:
            raise ValueError(f"not a fasta file: {str(b)}")
        e = m.start()
        desc = b[0:e]
        sid, _ = desc.split(b" ", maxsplit=1)
        seq = b[e + 1 :]
        seq = remove_white(seq)
        encoding = self.encoding
        return Record(
            sid.decode(encoding),
            desc.strip().decode(encoding),
            seq.upper().decode(encoding),
        )

    def __len__(self) -> int:
        return len(self.pos)

    @overload
    def __getitem__(self, idx: int) -> Record: ...
    @overload
    def __getitem__(self, idx: slice) -> list[Record]: ...
    @overload
    def __getitem__(self, idx: list[int]) -> list[Record]: ...
    def __getitem__(self, idx: int | slice | list[int]) -> Record | list[Record]:
        if isinstance(idx, int):
            return self.get_idx(idx)
        if isinstance(idx, list):
            return [self.get_idx(i) for i in idx]
        return [
            self.get_idx(i)
            for i in range(idx.start, idx.stop or len(self), idx.step or 1)
        ]


class CollectionFasta(Sequence[Record]):
    """Multiple memory mapped fasta files"""

    def __init__(self, fasta_files: Sequence[os.PathLike], encoding: str | None = None):
        self.fastas = [RandomFasta(f, encoding=encoding) for f in fasta_files]
        _cumsum = []
        cumsum = 0
        for f in self.fastas:
            cumsum += len(f)
            _cumsum.append(cumsum)
        self._cumsum = _cumsum
        self._len = cumsum

    def _map_idx(self, idx: int) -> tuple[int, RandomFasta]:
        i = bisect.bisect_right(self._cumsum, idx)
        if i > len(self._cumsum):
            raise IndexError("list out of range")
        r = self._cumsum[i - 1] if i > 0 else 0
        return idx - r, self.fastas[i]

    def _map_idxs(self, idxs: Sequence[int]) -> Iterator[tuple[int, RandomFasta]]:
        for idx in idxs:
            yield self._map_idx(idx)

    def __len__(self) -> int:
        return self._len

    def get_idx(self, idx: int) -> Record:
        i, rf = self._map_idx(idx)
        return rf.get_idx(i)

    def get_idxs(self, idxs: Sequence[int]) -> Iterator[Record]:
        for i, rf in self._map_idxs(idxs):
            yield rf.get_idx(i)

    @overload
    def __getitem__(self, idx: int) -> Record: ...
    @overload
    def __getitem__(self, idx: slice) -> list[Record]: ...
    @overload
    def __getitem__(self, idx: list[int]) -> list[Record]: ...
    def __getitem__(self, idx: int | slice | list[int]) -> Record | list[Record]:  # type: ignore
        if isinstance(idx, int):
            return self.get_idx(idx)
        if isinstance(idx, list):
            return [self.get_idx(i) for i in idx]
        return list(
            self.get_idxs(range(idx.start, idx.stop or len(self), idx.step or 1))
        )
