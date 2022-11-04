<p align="center">
 <img src="docs/images/logo-symbol-color.svg?raw=true", width="300", alt="SeqLike">
</p>

<h1 align="center"> SeqLike - flexible biological sequence objects in Python </h1>

<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
<p align="center">
   <a href="#contributors-">
    <img src="https://img.shields.io/badge/all_contributors-11-orange.svg?style=flat-square">
  </a>
  <a href="https://github.com/modernatx/seqlike/blob/main/LICENSE">
    <img src="https://img.shields.io/badge/license-Apache%202-blue">
  </a>
</p>
<!-- ALL-CONTRIBUTORS-BADGE:END -->

<!--- BADGES: START --->
<p align="center">
  <!--- This badge will work once supported python version information is available on PyPI --->
  <a href="https://pypi.org/project/seqlike/">
    <img src="https://img.shields.io/pypi/pyversions/seqlike?logo=pypi&style=flat&color=blue" alt="PyPI - Supported Python Version" />
  </a>
  <!--- Supported Python Version Badge End --->
  <a href="https://pypi.org/project/seqlike/">
    <img src="https://img.shields.io/pypi/v/seqlike?logo=pypi&style=flat&color=orange" alt="PyPI - Package Version" />
  </a>
  <a href="https://anaconda.org/conda-forge/seqlike">
    <img src="https://img.shields.io/conda/pn/conda-forge/seqlike?logo=anaconda&style=flat" alt="Conda - Platform" />
  </a>
  <a href="https://anaconda.org/conda-forge/seqlike">
    <img src="https://img.shields.io/conda/vn/conda-forge/seqlike?logo=anaconda&style=flat&color=orange" alt="Conda (channel only)" />
  </a>
  <a href="https://modernatx.github.io/seqlike">
    <img src="https://img.shields.io/static/v1?logo=github&style=flat&color=pink&label=docs&message=seqlike" alt="Docs - GitHub.io" />
  </a>
</p>
<!--- BADGES: END --->

## Introduction

A single object API that makes working with biological sequences in Python
 more ergonomic. It'll handle anything _like a sequence_.  

Built around the [Biopython SeqRecord class](https://biopython.org/wiki/SeqRecord),
SeqLikes abstract over the semantics of molecular biology (DNA -> RNA -> AA)
and data structures (strings, Seqs, SeqRecords, numerical encodings)
to allow manipulation of a biological sequence
at the level which is most computationally convenient.

## Code samples and examples

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

### Easily plot multiple sequence alignments

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

### Flexibly build and parse numerical sequence representations

```python
# Assume you have a dataframe with a column of 10 SeqLikes of length 90
df["seqs"].seq.to_onehot().shape # (10, 90, 23), padded if needed
```

To see more in action,
please check out the [docs](https://modernatx.github.io/seqlike/)!

<!-- ![Logo](https://dev-to-uploads.s3.amazonaws.com/uploads/articles/th5xamgrr6se0x5ro4g6.png) -->


## Getting Started

Install the library with `pip` or `conda`.

**With pip**

```python
pip install seqlike
```

**With conda**

```sh
conda install -c conda-forge seqlike
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

## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center"><a href="https://github.com/ndousis"><img src="https://avatars.githubusercontent.com/u/15198691?v=4?s=100" width="100px;" alt="Nasos Dousis"/><br /><sub><b>Nasos Dousis</b></sub></a><br /><a href="https://github.com/modernatx/seqlike/commits?author=ndousis" title="Code">💻</a></td>
      <td align="center"><a href="http://giessel.com"><img src="https://avatars.githubusercontent.com/u/1160997?v=4?s=100" width="100px;" alt="andrew giessel"/><br /><sub><b>andrew giessel</b></sub></a><br /><a href="https://github.com/modernatx/seqlike/commits?author=andrewgiessel" title="Code">💻</a></td>
      <td align="center"><a href="https://github.com/maxasauruswall"><img src="https://avatars.githubusercontent.com/u/14082213?v=4?s=100" width="100px;" alt="Max Wall"/><br /><sub><b>Max Wall</b></sub></a><br /><a href="https://github.com/modernatx/seqlike/commits?author=maxasauruswall" title="Code">💻</a> <a href="https://github.com/modernatx/seqlike/commits?author=maxasauruswall" title="Documentation">📖</a></td>
      <td align="center"><a href="https://ericmjl.github.io/"><img src="https://avatars.githubusercontent.com/u/2631566?v=4?s=100" width="100px;" alt="Eric Ma"/><br /><sub><b>Eric Ma</b></sub></a><br /><a href="https://github.com/modernatx/seqlike/commits?author=ericmjl" title="Code">💻</a> <a href="https://github.com/modernatx/seqlike/commits?author=ericmjl" title="Documentation">📖</a></td>
      <td align="center"><a href="https://github.com/MihirMetkar"><img src="https://avatars.githubusercontent.com/u/9938754?v=4?s=100" width="100px;" alt="Mihir Metkar"/><br /><sub><b>Mihir Metkar</b></sub></a><br /><a href="#ideas-MihirMetkar" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/modernatx/seqlike/commits?author=MihirMetkar" title="Code">💻</a></td>
      <td align="center"><a href="https://github.com/mccaron707"><img src="https://avatars.githubusercontent.com/u/26267127?v=4?s=100" width="100px;" alt="Marcus Caron"/><br /><sub><b>Marcus Caron</b></sub></a><br /><a href="https://github.com/modernatx/seqlike/commits?author=mccaron707" title="Documentation">📖</a></td>
      <td align="center"><a href="https://github.com/pagpires"><img src="https://avatars.githubusercontent.com/u/7856031?v=4?s=100" width="100px;" alt="pagpires"/><br /><sub><b>pagpires</b></sub></a><br /><a href="https://github.com/modernatx/seqlike/commits?author=pagpires" title="Documentation">📖</a></td>
    </tr>
    <tr>
      <td align="center"><a href="https://www.linkedin.com/in/sugatoray/"><img src="https://avatars.githubusercontent.com/u/10201242?v=4?s=100" width="100px;" alt="Sugato Ray"/><br /><sub><b>Sugato Ray</b></sub></a><br /><a href="#infra-sugatoray" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a> <a href="#maintenance-sugatoray" title="Maintenance">🚧</a></td>
      <td align="center"><a href="http://dmnfarrell.github.io/"><img src="https://avatars.githubusercontent.com/u/7859189?v=4?s=100" width="100px;" alt="Damien Farrell"/><br /><sub><b>Damien Farrell</b></sub></a><br /><a href="https://github.com/modernatx/seqlike/commits?author=dmnfarrell" title="Code">💻</a></td>
      <td align="center"><a href="https://github.com/farbod-nobar"><img src="https://avatars.githubusercontent.com/u/44842525?v=4?s=100" width="100px;" alt="Farbod Mahmoudinobar"/><br /><sub><b>Farbod Mahmoudinobar</b></sub></a><br /><a href="https://github.com/modernatx/seqlike/commits?author=farbod-nobar" title="Code">💻</a></td>
      <td align="center"><a href="https://github.com/JacobHayes"><img src="https://avatars.githubusercontent.com/u/2555532?v=4?s=100" width="100px;" alt="Jacob Hayes"/><br /><sub><b>Jacob Hayes</b></sub></a><br /><a href="#infra-JacobHayes" title="Infrastructure (Hosting, Build-Tools, etc)">🚇</a></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!
