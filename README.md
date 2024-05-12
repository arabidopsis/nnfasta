# nnfasta

Lightweight Neural Net efficient FASTA dataset suitable for training.

Should be memory efficient across process boundaries.
So useful as input to torch/tensorflow dataloaders with multiple workers etc.
(see [this issue](https://github.com/pytorch/pytorch/issues/13246#issuecomment-905703662))

Presents a list of fasta files as a simple `abc.Sequence`
so you can inquire about `len(dataset)` and retrieve
`Record`s randomly with `dataset[i]`

Uses Python's `mmap.mmap` under the hood.

The underlying FASTA's should be "well formed" since there is
minimal sanity checking done.

## Install

Install with:

```bash
pip install nnfasta
```

**There are no dependencies**, you just need a modern (>= 3.9) python (< 12K of code).

## Usage

```python
from nnfasta import nnfastas

dataset = nnfastas(['athaliana.fasta','triticum.fasta','zmays.fasta'])

# display the number of sequences
print(len(dataset))

# get a particular record
rec = dataset[20]
print('sequence', rec.id, rec.description, rec.seq)
```

**Warning**: No checks are made for the existence of
the fasta files. Also files of zero length will be rejected
by `mmap`.

A `Record` mimics biopython's [`SeqRecord`](https://biopython.org/wiki/SeqRecord) and is simply:

```python
@dataclass
class Record:
    id: str
    """Sequence ID"""
    description: str
    """Line prefixed by '>'"""
    seq: str
    """Sequence stripped of whitespace and uppercased"""

    @property
    def name(self) -> str:
        return self.id
```

The major difference is that `seq` is just a simple `str` not a biopython `Seq` object
(We just don't want the `Bio` dependency -- `nnfasta` has _no_ dependencies).

## Arguments

You can give `nnfastas` either a filename, a `Path`, the actual
bytes in the file or an open file pointer (opened with `mode="rb"`)
_OR_ a list of these things. e.g:

```python

from nnfasta import nnfastas
my = "my.fasta"
fa = nnfastas([my, open(my mode="rb"),
            open(my, mode="rb").read()])
```

## Encoding

The files are assumed to be encoded as "`ASCII`". If this is not the
case then `nnfastas` accepts an `encoding` argument. All the files
presented to `nnfastas` are assumed to be similarly encoded. You can
alter the decoding with the `errors` keyword (default=`strict`).

## Test and Train Split best practice

Use `SubsetFasta`

```python
from nnfasta import nnfastas, SubsetFasta
from sklearn.model_selection import train_test_split

dataset = nnfastas(['athaliana.fasta','triticum.fasta','zmays.fasta'])
train_idx, test_idx = train_test_split(range(len(dataset)),test_size=.1,shuffle=True)

# these are still Sequence[Record] objects.

train_data = SubsetFasta(dataset, train_idx)
test_data = SubsetFasta(dataset, test_idx)

# *OR* ... this is basically the same
import torch
train_data, test_data = torch.utils.data.random_split(dataset, [.9, .1])

```

See the pytorch `Subset` logic [here](https://pytorch.org/docs/stable/data.html#torch.utils.data.Subset)


## How it works

We memory map the input files and use python's `re` package to scan the files
for `b"\r>|\n>|^>"`  bytes from which we compute a start, end for each
record and create an `array.array` (in memory).

The operating system will ensure that similar mmapped pages in different process
are shared.

Enjoy peps!
