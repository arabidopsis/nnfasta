import bisect
import mmap
import re
import os
import io
from collections.abc import Sequence
from dataclasses import dataclass, replace
from typing import Iterator, overload, TypeAlias, cast, IO, Any
import array

# a "fasta" file is a Path or filename or the actual bytes or an
# open file
Fasta: TypeAlias = os.PathLike | str | bytes | IO[bytes]

# This mimics biopython's SeqRecord. Except that `seq`
# is just a string and not a `Seq`
@dataclass
class Record:
    """Mimics biopython Record"""

    id: str
    """Sequence ID"""
    description: str
    """Line prefixed by '>'"""
    seq: str
    """Sequence stripped of whitespace and uppercased"""

    # just fake SeqRecord
    # still missing count(sub) translate() reverse_complement()
    @property
    def name(self) -> str:
        """same as ID"""
        return self.id

    @property
    def features(self) -> list[Any]:
        """SeqRecord list of features"""
        return []

    @property
    def dbxrefs(self) -> list[Any]:
        """SeqRecord list of dbxrefs"""
        return []

    @property
    def annotations(self) -> dict[str, Any]:
        """SeqRecord dict of annotations"""
        return {}

    @property
    def letter_annotations(self) -> dict[str, Any]:
        """SeqRecord letter_annotations"""
        return {}

    def format(self, fmt: str = "fasta") -> str:
        """format record as FASTA"""
        from itertools import batched

        if fmt.lower() != "fasta":
            raise ValueError("can only format as FASTA")
        s = "\n".join("".join(s) for s in batched(self.seq, 80))
        return f">{self.description}\n{s}"

    def islower(self) -> bool:
        """Is sequence in lowercase"""
        return self.seq.lower() == self.seq

    def isupper(self) -> bool:
        """is sequence in uppercalse"""
        return self.seq.upper() == self.seq

    def upper(self) -> "Record":
        """Uppercase sequence"""
        return replace(self, seq=self.seq.upper())

    def lower(self) -> "Record":
        """Lowercase sequence"""
        return replace(self, seq=self.seq.lower())


PREFIX = re.compile(b"(\n>|\r>|^>)", re.M)

WHITE = re.compile(rb"\W+", re.I | re.M)
EOL = re.compile(rb"\n|\r", re.I | re.M)


def remove_white(s: bytes) -> bytes:
    """remove whitespace from byte list"""
    return WHITE.sub(b"", s)


def isopen(obj: Sequence[Fasta] | Fasta) -> bool:
    """check if an object is an open file or not"""
    # return hasattr(obj, 'close')
    return isinstance(obj, io.IOBase)


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
        text encoding of these files [default: ascii]
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

    if isinstance(fasta_file_or_bytes, (os.PathLike, str, bytes, io.IOBase)):
        fasta_file_or_bytes = [fasta_file_or_bytes]  # type: ignore
    assert isinstance(fasta_file_or_bytes, Sequence)
    if len(fasta_file_or_bytes) == 1:
        return RandomFasta(fasta_file_or_bytes[0], encoding=encoding, errors=errors)  # type: ignore
    return CollectionFasta(fasta_file_or_bytes, encoding=encoding, errors=errors)  # type: ignore


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
        if isinstance(fasta_file_or_bytes, bytes):
            self.isopen = False
            self.fasta = fasta_file_or_bytes
            self.fp = None
        else:
            self.isopen = isopen(fasta_file_or_bytes)
            if self.isopen:
                self.fp = cast(IO[bytes], fasta_file_or_bytes)
            else:
                self.fp = open(fasta_file_or_bytes, "rb")  # type: ignore
            self.fasta = cast(
                bytes, mmap.mmap(self.fp.fileno(), 0, prot=mmap.PROT_READ)
            )
        self._pos = self._find_pos()

    def __del__(self):
        if self is not None:
            if self.fp and not self.isopen:
                self.fp.close()
            self.fp = None

    def _find_pos(self) -> array.ArrayType:
        end, start = zip(*((h.start(), h.end()) for h in PREFIX.finditer(self.fasta)))

        end = end[1:] + (len(self.fasta),)
        # we return an array so that
        # we have fewer refcounts
        return array.array("Q", [a for se in zip(start, end) for a in se])

    def get_idx(self, idx: int) -> Record:
        """get Record for index"""
        if idx < 0:
            idx = len(self) - ((-idx) % len(self))
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
        return [self.get_idx(i) for i in idx]

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
        assert len(self.fastas) > 0, "list of fasta files should not be empty"
        _cumsum = []
        cumsum = 0
        for f in self.fastas:
            cumsum += len(f)
            _cumsum.append(cumsum)
        self._cumsum = _cumsum

    def _map_idx(self, idx: int) -> tuple[int, RandomFasta]:
        if idx < 0:
            if -idx > len(self):
                raise ValueError(
                    "absolute value of index should not exceed dataset length"
                )
            idx = len(self) + idx
        i = bisect.bisect_right(self._cumsum, idx)
        if i > len(self._cumsum):
            raise IndexError("list out of range")
        r = self._cumsum[i - 1] if i > 0 else 0
        return idx - r, self.fastas[i]

    def _map_idxs(self, idxs: Sequence[int]) -> Iterator[tuple[int, RandomFasta]]:
        for idx in idxs:
            yield self._map_idx(idx)

    def __len__(self) -> int:
        return self._cumsum[-1]

    def get_idx(self, idx: int) -> Record:
        """Given an integer ID return the Record"""
        i, rf = self._map_idx(idx)
        return rf.get_idx(i)

    def get_idxs(self, idxs: Sequence[int]) -> Iterator[Record]:
        """get Records for sequence of integers"""
        for i, rf in self._map_idxs(idxs):
            yield rf.get_idx(i)

    def __getitems__(self, idx: list[int]) -> list[Record]:
        """torch extention"""
        return [self.get_idx(i) for i in idx]

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


class SubsetFasta(Sequence[Record]):
    """Return sequence records from a list of indicies.

    See https://pytorch.org/docs/stable/data.html#torch.utils.data.Subset
    """

    def __init__(self, dataset: Sequence[Record], indexes: Sequence[int]):
        """Use index to create a new dataset from another"""
        self._dataset = dataset
        self._indexes = indexes
        assert len(dataset) > max(self._indexes) and min(self._indexes) >= 0

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
    def __getitem__(self, idx: int | slice | list[int]) -> Record | list[Record]:  # type: ignore
        index = self._indexes
        if isinstance(idx, int):
            return self._dataset[index[idx]]
        if isinstance(idx, list):
            return [self._dataset[index[i]] for i in idx]
        return [
            self._dataset[index[i]]
            for i in range(idx.start, idx.stop or len(self), idx.step or 1)
        ]
