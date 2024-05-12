from __future__ import annotations

import array
import bisect
import io
import mmap
import os
import re
from collections.abc import Sequence
from dataclasses import dataclass
from itertools import islice
from typing import cast
from typing import Iterable
from typing import Iterator
from typing import overload
from typing import TypeAlias
from typing import TypeVar

T = TypeVar("T")


# batched backport. Defined in itertools in 3.12
def batched(iterable: Iterable[T], n: int) -> Iterator[tuple[T, ...]]:
    """Batch incoming iterable"""
    # batched('ABCDEFG', 3) → ABC DEF G
    if n < 1:
        raise ValueError("n must be at least one")
    it = iter(iterable)
    while True:
        batch = tuple(islice(it, n))
        if not batch:
            break
        yield batch


# a "fasta" file is a Path or filename or the actual bytes or an
# open file
# WARNING isinstance("string", Sequence) is True
# and isinstance(b"string", Sequence) is True
Fasta: TypeAlias = os.PathLike | str | bytes | io.IOBase
FType = (os.PathLike, str, bytes, io.IOBase)


# This mimics biopython's SeqRecord. Except that `seq`
# is just a string and not a `Seq`
@dataclass
class Record:
    """Mimics biopython SeqRecord.

    The main difference is that `seq` is just a simple string.
    The class is also missing `translate()` and `reverse_complement()`
    functions. Really get biopython if you need these!

    ``from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    rec = SeqRecord(id=seq.id, description=rec.description, seq=Seq(rec.seq))``
    """

    id: str
    """Sequence ID"""
    description: str
    """Line prefixed by '>'"""
    seq: str
    """Sequence stripped of whitespace and uppercased"""

    @property
    def name(self) -> str:
        """same as ID"""
        return self.id

    # just fake SeqRecord
    # still missing translate() reverse_complement()
    # and dunder calls
    # def to_rec(self):
    #     """Return a biopython SeqRecord (if biopython is installed)"""
    #     from Bio.SeqRecord import SeqRecord
    #     from Bio.Seq import Seq

    #     return SeqRecord(id=self.id, description=self.description, seq=Seq(self.seq))

    # @property
    # def features(self) -> list[Any]:
    #     """SeqRecord list of features"""
    #     return []

    # @property
    # def dbxrefs(self) -> list[Any]:
    #     """SeqRecord list of dbxrefs"""
    #     return []

    # @property
    # def annotations(self) -> dict[str, Any]:
    #     """SeqRecord dict of annotations"""
    #     return {}

    # @property
    # def letter_annotations(self) -> dict[str, Any]:
    #     """SeqRecord letter_annotations"""
    #     return {}

    def format(self, fmt: str = "fasta") -> str:
        """format record as FASTA"""
        return self.__format__(fmt)  # pylint: disable=unnecessary-dunder-call

    def __format__(self, fmt: str) -> str:
        """format record as FASTA"""

        if fmt.lower() != "fasta":
            raise ValueError("can only format as FASTA")
        s = "\n".join("".join(s) for s in batched(self.seq, 60))
        return f">{self.description}\n{s}\n"

    # def islower(self) -> bool:
    #     """Is sequence in lowercase"""
    #     return self.seq.lower() == self.seq

    # def isupper(self) -> bool:
    #     """is sequence in uppercalse"""
    #     return self.seq.upper() == self.seq

    # def upper(self) -> Record:
    #     """Uppercase sequence"""
    #     return replace(self, seq=self.seq.upper())

    # def lower(self) -> Record:
    #     """Lowercase sequence"""
    #     return replace(self, seq=self.seq.lower())

    # def count(self, sub: bytes | str, start: int | None = None, end: int | None = None):
    #     """Return the number of non-overlapping occurrences of sub in data[start:end].

    #     Optional arguments start and end are interpreted as in slice notation.
    #     This method behaves as the count method of Python strings.
    #     """
    #     if isinstance(sub, str):
    #         sub = sub.encode("ASCII")
    #     return self.seq.encode("ASCII").count(sub, start, end)


PREFIX = re.compile(b"(\n>|\r>|^>)", re.M)

WHITE = re.compile(rb"\W+", re.I | re.M)
EOL = re.compile(rb"\n|\r", re.I | re.M)


def remove_white(s: bytes) -> bytes:
    """remove whitespace from byte list"""
    return WHITE.sub(b"", s)


def nnfastas(
    fasta_file_or_bytes: Sequence[Fasta] | Fasta,
    *,
    encoding: str | None = None,
    errors: str | None = None,
) -> Sequence[Record]:
    """Given a sequence of fasta files return an indexable (list like) Fasta object.

    Parameters
    ----------

    fasta_files: Sequence[PathLike | str | bytes | IO[bytes]] | (PathLike| str | bytes | IO[bytes])
        sequence of Fasta files to mmap. (Can be the raw bytes from the file too or the 'rb' opened file!)
    encoding: str, optional
        text encoding of these files [default: 'ascii']
    errors: str, optional
        how to treat decoding errors. default='strict'

    Returns
    -------

    A ``Sequence[Record]`` object.
    """
    # WARNING isinstance("string", Sequence) is True
    # and isinstance(b"string", Sequence) is True
    if not fasta_file_or_bytes:
        raise ValueError("no fasta files!")

    # current mypy just does not get it...
    if isinstance(fasta_file_or_bytes, FType):
        fasta_file_or_bytes = [fasta_file_or_bytes]
    assert isinstance(fasta_file_or_bytes, Sequence) and not isinstance(
        fasta_file_or_bytes,
        bytes,
    )
    if len(fasta_file_or_bytes) == 1:
        return RandomFasta(fasta_file_or_bytes[0], encoding=encoding, errors=errors)
    return CollectionFasta(fasta_file_or_bytes, encoding=encoding, errors=errors)


class RandomFasta(Sequence[Record]):
    """Memory mapped fasta file."""

    ENCODING = "ascii"
    ERRORS = "strict"

    def __init__(
        self,
        fasta_file_or_bytes: Fasta,
        *,
        encoding: str | None = None,
        errors: str | None = None,
    ):

        self.encoding = encoding or self.ENCODING
        self.errors = errors or self.ERRORS
        self._fp: io.IOBase | None
        if isinstance(fasta_file_or_bytes, bytes):
            self.isopen = False
            self.fasta = fasta_file_or_bytes
            self._fp = None
        else:
            isopen = isinstance(fasta_file_or_bytes, io.IOBase)
            if isopen:
                assert isinstance(fasta_file_or_bytes, io.IOBase)  # mypy you dummy!
                self._fp = fasta_file_or_bytes
            else:
                assert not isinstance(fasta_file_or_bytes, (io.IOBase, bytes))
                self._fp = open(fasta_file_or_bytes, mode="rb")
            assert self._fp is not None
            self.isopen = isopen
            # we cast to bytes to ensure we don't use any mmap functions....
            self.fasta = cast(
                bytes,
                mmap.mmap(self._fp.fileno(), 0, access=mmap.ACCESS_READ),
            )
        self._pos = self._find_pos()

    def __del__(self):
        if self is not None:
            if self._fp and not self.isopen:
                self._fp.close()
            self._fp = None

    def _find_pos(self) -> array.ArrayType:
        end, start = zip(*((h.start(), h.end()) for h in PREFIX.finditer(self.fasta)))
        end = end[1:] + (len(self.fasta),)
        # we return an array so that
        # we have fewer refcounts
        return array.array("Q", [a for se in zip(start, end) for a in se])

    def _get_idx(self, idx: int) -> Record:
        """get Record for index"""
        if idx < 0:
            n = len(self)
            if idx < -n:
                raise IndexError("index out of range")
            idx = n + idx
            # idx = len(self) - ((-idx) % len(self))

        idx = 2 * idx
        s, e = self._pos[idx], self._pos[idx + 1]
        b = self.fasta[s:e]  # mmap goes to disk
        m = EOL.search(b)
        if not m:
            raise ValueError(f"not a fasta file: {str(b)}")
        e = m.start()
        desc = b[:e]
        if b" " in desc:
            sid, _ = desc.split(b" ", maxsplit=1)
        else:
            sid = desc
        seq = b[e + 1 :]
        seq = remove_white(seq)
        encoding = self.encoding
        errors = self.errors
        return Record(
            sid.decode(encoding, errors=errors),
            desc.strip().decode(encoding, errors=errors),
            seq.upper().decode(encoding, errors=errors),
        )

    def __len__(self) -> int:
        return len(self._pos) // 2

    def __getitems__(self, idx: list[int]) -> list[Record]:
        """torch extention"""
        return [self._get_idx(i) for i in idx]

    @overload
    def __getitem__(self, idx: int) -> Record: ...
    @overload
    def __getitem__(self, idx: slice) -> list[Record]: ...
    @overload
    def __getitem__(self, idx: list[int]) -> list[Record]: ...

    def __getitem__(self, idx: int | slice | list[int]) -> Record | list[Record]:
        if isinstance(idx, int):
            return self._get_idx(idx)
        if isinstance(idx, list):
            return [self._get_idx(i) for i in idx]
        return [
            self._get_idx(i)
            for i in range(idx.start, idx.stop or len(self), idx.step or 1)
        ]

    def close(self) -> None:
        """Close any open files"""
        if self and self._fp is not None:
            self._fp.close()
            self._fp = None


class CollectionFasta(Sequence[Record]):
    """Multiple memory mapped fasta files"""

    # see also https://pytorch.org/docs/stable/_modules/torch/utils/data/dataset.html#ConcatDataset

    def __init__(
        self,
        fasta_file_or_bytes: Sequence[Fasta],
        *,
        encoding: str | None = None,
        errors: str | None = None,
    ):
        self.fastas = [
            RandomFasta(f, encoding=encoding, errors=errors)
            for f in fasta_file_or_bytes
        ]
        if len(self.fastas) == 0:
            raise ValueError("list of fasta files should not be empty")
        _cumsum = []
        cumsum = 0
        for f in self.fastas:
            cumsum += len(f)
            _cumsum.append(cumsum)
        self._cumsum = array.array("Q", _cumsum)

    def _map_idx(self, idx: int) -> tuple[int, RandomFasta]:
        if idx < 0:
            n = len(self)
            if idx < -n:
                raise IndexError("index out of range")
            idx = n + idx
            # idx = len(self) - ((-idx) % len(self))
        i = bisect.bisect_right(self._cumsum, idx)
        if i > len(self._cumsum):
            raise IndexError("index out of range")
        r = self._cumsum[i - 1] if i > 0 else 0
        return idx - r, self.fastas[i]

    def _map_idxs(self, idxs: Sequence[int]) -> Iterator[tuple[int, RandomFasta]]:
        for idx in idxs:
            yield self._map_idx(idx)

    def __len__(self) -> int:
        return self._cumsum[-1]

    def _get_idx(self, idx: int) -> Record:
        """Given an integer ID return the Record"""
        i, rf = self._map_idx(idx)
        return rf._get_idx(i)  # pylint: disable=protected-access

    def _get_idxs(self, idxs: Sequence[int]) -> Iterator[Record]:
        """get Records for sequence of integers"""
        for i, rf in self._map_idxs(idxs):
            yield rf._get_idx(i)  # pylint: disable=protected-access

    def __getitems__(self, idx: list[int]) -> list[Record]:
        """torch extention"""
        return [self._get_idx(i) for i in idx]

    @overload
    def __getitem__(self, idx: int) -> Record: ...
    @overload
    def __getitem__(self, idx: slice) -> list[Record]: ...
    @overload
    def __getitem__(self, idx: list[int]) -> list[Record]: ...

    def __getitem__(self, idx: int | slice | list[int]) -> Record | list[Record]:
        if isinstance(idx, int):
            return self._get_idx(idx)
        if isinstance(idx, list):
            return [self._get_idx(i) for i in idx]
        return list(
            self._get_idxs(range(idx.start, idx.stop or len(self), idx.step or 1)),
        )

    def close(self) -> None:
        """Close any open files"""
        if self:
            for b in self.fastas:
                b.close()


class SubsetFasta(Sequence[Record]):
    """Return sequence records from a list of indicies.

    See https://pytorch.org/docs/stable/data.html#torch.utils.data.Subset
    """

    def __init__(self, dataset: Sequence[Record], indexes: Sequence[int]):
        """Use index to create a new dataset from another"""
        self._dataset = dataset
        self._indexes = indexes
        if len(dataset) <= max(self._indexes) or min(self._indexes) < 0:
            raise ValueError("indexes don't index dataset")

    def __len__(self) -> int:
        return len(self._indexes)

    def __getitems__(self, idx: list[int]) -> list[Record]:
        """torch extention"""
        index = self._indexes
        return [self._dataset[index[i]] for i in idx]

    @overload
    def __getitem__(self, idx: int) -> Record: ...
    @overload
    def __getitem__(self, idx: slice) -> list[Record]: ...
    @overload
    def __getitem__(self, idx: list[int]) -> list[Record]: ...

    def __getitem__(self, idx: int | slice | list[int]) -> Record | list[Record]:
        index = self._indexes
        if isinstance(idx, int):
            return self._dataset[index[idx]]
        if isinstance(idx, list):
            return [self._dataset[index[i]] for i in idx]
        return [
            self._dataset[index[i]]
            for i in range(idx.start, idx.stop or len(self), idx.step or 1)
        ]
