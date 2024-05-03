# nnfasta

Neural Net efficient fasta Dataset for Training.

Should be efficient across process boundaries etc.
So useful as input to torch/tensorflow dataloaders etc.

Presents a list of fasta files as a simple `abc.Sequence`
so you can inquire about `len(dataset)` and retrieve
`Record`s with `dataset[i]`

## Usage

Install:

```bash
pip install nnfasta
```

```python

from nnfasta import nnfastas 


dataset = nnfastas(['athaliana.fasta','triticum.fasta','zmays.fasta'])`

print(len(datset))

rec = dataset[20]
print('sequence', rec.seq)
```
