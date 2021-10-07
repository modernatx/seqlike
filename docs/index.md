# SeqLike: Flexible Biological Sequence Objects in Python

SeqLike is a Python package that lets you conveniently manipulate biological sequences.
It solves some of the following problems:

1. Sequence representation inter-conversion (AA vs. NT, and str/Seq/SeqRecord/arrays) via a single object's API.
2. Processing a collection of sequences easily in Python without needing to switch out to a shell.
3. Convenience APIs to visualize of a collection of sequences.
## Installation

You can install SeqLike from PyPI:

```bash
pip install seqlike
```

SeqLike can also be installed from conda-forge

```bash
conda install -c conda-forge seqlike
```

## Usage/Examples

Please see the 5-minute tutorial notebook for a quick guide on how to use SeqLike.
## Features

- Interconversion between AA and NT forms.
- Pandas accessor methods for manipulating collections of sequences.
- Convenient multiple sequence alignment and plotting without switching to a shell.
- BioPython SeqRecord behaviour.

## License

SeqLike is licensed under the [Apache 2.0 license](https://choosealicense.com/licenses/apache-2.0/).
