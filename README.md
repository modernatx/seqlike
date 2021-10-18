# SeqLike

[![Apache v2 License](https://img.shields.io/badge/license-Apache%202-blue)](https://github.com/modernatx/seqlike/blob/main/LICENSE)

A single object API that makes working with biological sequences in Python
just a bit more ergonomic.

We introduce the concept of a "Sequence" as a consistent abstract entity that is maintained accross biological (nucleotide -> protein) and computational (string -> encoding) transforms.

The SeqLike class includes many convenience methods
for converting between representations:
One-Hot vs Index Encodings,
strings,
Biopython Seqs,
etc.

Built around the [Biopython SeqRecord class](https://biopython.org/wiki/SeqRecord),
SeqLikes abstracts over the semantics of molecular biology (DNA -> RNA -> AA)
to allow manipulating an abstract sequence
at the level which is most computationally convenient.

For example:

- Perform multiple sequence alignment and visualize the result [without leaving your Jupyter Notebook](https://github.com/modernatx/seqlike/blob/main/notebooks/tutorial.ipynb).
- Move back and forth between nucleotide and amino-acid sequences while preserving their biological semantics.

To see more in action,
please check out the [docs](https://modernatx.github.io/seqlike/)!

<!-- ![Logo](https://dev-to-uploads.s3.amazonaws.com/uploads/articles/th5xamgrr6se0x5ro4g6.png) -->

## Code samples

### Build data-type agnostic functions

```python
def f(seq: SeqLikeType, *args):
	seq = SeqLike(seq, seq_type="nt").to_seqrecord()
	# ...
```

#### Streamline conversion to/from ML friendly representations

```python
prediction = model(aaSeqLike('MSKGEELFTG').to_onehot())
new_seq = ntSeqLike(generative_model.sample(), alphabet="-ACGTUN")
```

### Interconvert between AA and NT forms of a sequence

Back-translation is conveniently built-in!

```python
s_nt = ntSeqLike("ATGTCTAAAGGTGAA")
s_nt[0:3] # ATG
s_nt.aa()[0:3] # MSK, nt->aa is well defined
s_nt.aa()[0:3].nt() # ATGTCTAAA, works because SeqLike now has both reps
s_nt[:-1].aa() # TypeError, len(s_nt) not a multiple of 3

s_aa = aaSeqLike("MSKGE")
s_aa.nt() # AttributeError, aa->nt is undefined w/o codon map
s_aa = aaSeqLike(s_aa, codon_map=random_codon_map)
s_aa.nt() # now works, backtranslated to e.g. ATGTCTAAAGGTGAA
s_aa[:1].nt() # ATG, codon_map is maintained
```

### Plot the multiple sequence alignment

```python
seqs = [s for s in SeqIO.parse("file.fasta", "fasta")]
df = pd.DataFrame(
    {
        "names": [s.name for s in seqs],
        "seqs": [aaSeqLike(s) for s in seqs],
    }
)
df["aligned"] = df["seqs"].seq.align()
df["aligned"].seq.plot()
```

### Easily build numerical representations

```python
# Assume you have a dataframe with a column of SeqLikes
df["seqs"].seq.to_onehot().shape # (10, 91, 23), padded if needed
```

## Getting Started

```python
pip install seqlike
```

## Authors

- [@andrewgiessel](https://github.com/andrewgiessel)
- [@maxasauruswall](https://github.com/maxasauruswall)
- [@MihirMetkar](https://github.com/MihirMetkar)
- [@ndousis](https://github.com/ndousis)
- [@ericmjl](https://github.com/ericmjl)

## Support

- Questions about usage should be posed on [Stack Overflow with the #seqlike tag][SO].
- Bug reports and feature requests are managed using the [Github issue tracker][gh_issues].

[SO]: https://stackoverflow.com/questions/tagged/seqlike
[gh_issues]: https://github.com/modernatx/seqlike/issues
