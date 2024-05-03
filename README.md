# nnfasta

Neural Net efficient fasta Dataset for Training.

Should be memory efficient across process boundaries.
So useful as input to torch/tensorflow dataloaders etc.

Presents a list of fasta files as a simple `abc.Sequence`
so you can inquire about `len(dataset)` and retrieve
`Record`s with `dataset[i]`

## Install

Install:

```bash
pip install nnfasta
```

There are **no** dependencies.

## Usage

```python

from nnfasta import nnfastas 


dataset = nnfastas(['athaliana.fasta','triticum.fasta','zmays.fasta'])

# display number of sequences
print(len(dataset))

# get a particular record
rec = dataset[20]
print('sequence', rec.id, rec.description, rec.seq)
```
