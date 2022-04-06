from copy import deepcopy
import tempfile


import numpy as np
import pandas as pd
import pytest
from PIL import Image

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from seqlike import SeqLike
from seqlike.codon_tables import human_codon_table, human_codon_map, codon_table_to_codon_map

from . import test_path


# TODO: Turn this into a pytest fixture using Hypothesis.
# We might need to refactor out the fixtures a bit.
nt_seqs = [SeqLike(s, "nt") for s in SeqIO.parse(test_path / f"abs_nt_4.fasta", "fasta")]
s = SeqLike(SeqIO.read(test_path / f"test.fa", "fasta"), seq_type="dna")
s_aa = SeqLike(SeqIO.read(test_path / f"test.fa", "fasta"), seq_type="dna").aa()
s_aa_with_codon_map = SeqLike(
    SeqIO.read(test_path / f"test.fa", "fasta"),
    codon_map=human_codon_map,
    seq_type="dna",
).aa()

seqs = [deepcopy(s)] * 10
seqs_aa = [deepcopy(s_aa)] * 10
seqs_aa_with_codon_map = [deepcopy(s_aa_with_codon_map)] * 10
seqs_mixed = deepcopy(seqs) + deepcopy(seqs_aa)

# ---- test list of seqs of various types ---------
seqs_list = [
    (seqs, "nt"),
    (seqs_aa, "aa"),
    (seqs_aa_with_codon_map, "aa"),
    pytest.param(
        seqs_mixed,
        None,
        marks=pytest.mark.xfail(reason="Not a homogeneous list of SeqLikes."),
    ),
]


@pytest.mark.parametrize("seqs, _type", seqs_list)
def test_init_and__type(seqs, _type):
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": seqs})
    assert isinstance(df, pd.DataFrame)
    assert df.seqs.seq._type == _type.upper()


@pytest.mark.parametrize("seqs, _type", seqs_list)
def test_write(seqs, _type):
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": seqs})
    # r+ is read-writable
    with tempfile.NamedTemporaryFile(mode="r+") as tempf:
        df["seqs"].seq.write(tempf, "fasta")
        # rewind file after writing
        tempf.seek(0)
        read_seqs = pd.Series(SeqLike(s, seq_type=_type) for s in SeqIO.parse(tempf, "fasta"))
        for seq1, seq2 in zip(read_seqs, df["seqs"]):
            assert str(seq1) == str(seq2)


@pytest.mark.parametrize("seqs, _type", seqs_list)
def test_plot(seqs, _type):
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": seqs})
    assert isinstance(df.seqs.seq.plot(use_bokeh=False), Image.Image)


@pytest.mark.parametrize("seqs, _type", seqs_list)
def test_align(seqs, _type):
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": seqs})
    assert isinstance(df.seqs.seq.align(), pd.Series)


@pytest.mark.parametrize("seqs, _type", seqs_list)
def test_as_alignment(seqs, _type):
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": seqs})
    assert isinstance(df.seqs.seq.as_alignment(), MultipleSeqAlignment)


@pytest.mark.parametrize("seqs, _type", seqs_list)
def test_as_counts(seqs, _type):
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": seqs})
    as_counts = df.seqs.seq.as_counts()
    assert isinstance(as_counts, np.ndarray)
    assert as_counts.shape == (max(len(s) for s in seqs), len(seqs[0].alphabet))


@pytest.mark.parametrize("seqs, _type", seqs_list)
def test_extend_ambiguous_counts(seqs, _type):
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": seqs})
    extended_counts = df.seqs.seq._extend_ambiguous_counts()
    assert isinstance(extended_counts, np.ndarray)


@pytest.mark.parametrize("seqs, _type", seqs_list)
def test_consensus(seqs, _type):
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": seqs})
    consensus = df.seqs.seq.consensus()
    assert isinstance(consensus, SeqLike)
    assert len(consensus) == max(len(s) for s in seqs)


@pytest.mark.parametrize("seqs, _type", seqs_list)
def test_degenerate(seqs, _type):
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": seqs})
    degenerate = df.seqs.seq.degenerate()
    assert isinstance(degenerate, SeqLike)
    assert len(degenerate) == max(len(s) for s in seqs)
    assert set(degenerate).issubset(set(df.seqs.seq.alphabet))


def test_consensus2():
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": nt_seqs})
    consensus = df.seqs.seq.consensus()
    assert (
        str(consensus)
        == "TCAATTGGGGGAGGAGCTCTGGTGGAGGCGGTAGCGGAGGCGGAGGGTCGGCTAGCCAAGTCCAATTGGTTGAATCTGGTGGTGGTGTTGTTCAACCAGGTGGTTCTTTGAGATTGTCTT"
    )


def test_degenerate2():
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": nt_seqs})
    degenerate = df.seqs.seq.degenerate()
    assert (
        str(degenerate)
        == "TCAATTGGGGGAGGAGCTCTSGTGGWGGCVGTAGCGGAGKCGGAGGKTCSGCWAGCCAAGTCCAATTGGTTGAATCTGGTGGTGGTGTTGTTCAACCAGGTGGTTCTTTGAGATTGTCTT"
    )


@pytest.mark.parametrize("seqs, _type", seqs_list)
def test_max_length(seqs, _type):
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": seqs})
    assert df.seqs.seq.max_length() == max(len(x) for x in seqs)


@pytest.mark.parametrize("seqs, _type", seqs_list)
def test_get_seq_by_id(seqs, _type):
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": seqs})
    assert df.seqs.seq.get_seq_by_id(seqs[0].id) == seqs[0]


@pytest.mark.parametrize("seqs, _type", seqs_list)
def test__getitem__(seqs, _type):
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": seqs})
    assert isinstance(df.seqs.seq[:, :5], pd.Series)
    assert len(df.seqs.seq[:, :5]) == len(df)
    assert all([len(s) == 5 for s in df.seqs.seq[:, :5]])

    assert isinstance(df.seqs.seq[:2, :], pd.Series)
    assert len(df.seqs.seq[:2, :]) == 2

    assert isinstance(df.seqs.seq[0, :], SeqLike)
    assert isinstance(df.seqs.seq[0:1, :], pd.Series)
    assert len(df.seqs.seq[0:1, :]) == 1
    assert isinstance(df.seqs.seq[:, 0], pd.Series)
    assert len(df.seqs.seq[:, 0]) == len(df)

    assert df.seqs.seq.alphabet == df.seqs[1:].seq.alphabet


@pytest.mark.parametrize("seqs, _type", seqs_list)
def test_nt(seqs, _type):
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": seqs})
    assert isinstance(df.seqs.seq.nt(), pd.Series)
    assert all([s._type == "NT" for s in df.seqs.seq.nt()])


@pytest.mark.parametrize("seqs, _type", seqs_list)
def test_aa(seqs, _type):
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": seqs})
    assert isinstance(df.seqs.seq.aa(), pd.Series)
    assert all([s._type == "AA" for s in df.seqs.seq.aa()])


@pytest.mark.parametrize("seqs, _type", seqs_list)
def test_backtranslate(seqs, _type):
    # TODO: Refactor this test,
    # because it's covering too many different types of cases.

    df = pd.DataFrame({"seqs": seqs})

    # NOTE: I believe the test here is incorrect.
    # We should be testing that an AttributeError is raised.
    #
    # TODO: Check with Andrew.
    if _type == "aa" and df.seqs.iloc[0].codon_map is None:
        # The old implementation.
        # CAN't backtrans if no codon_map anywhere
        # assert isinstance(df.seqs.seq.back_translate(), pd.Series)

        # The proposed implementation
        with pytest.raises(AttributeError):
            df.seqs.seq.back_translate()

    elif _type == "aa" and df.seqs.iloc[0].codon_map is None:  # CAN backtrans if codon_map specified
        assert isinstance(
            df.seqs.seq.back_translate(codon_map=codon_table_to_codon_map(human_codon_table)),
            pd.Series,
        )

    elif (
        _type == "aa" and df.seqs.iloc[0].codon_map is not None
    ):  # CAN backtrans if codon_map is in the SeqLikes already
        assert isinstance(df.seqs.seq.back_translate(), pd.Series)

    elif _type == "nt":  # can't backtranslate from NT
        with pytest.raises(TypeError):
            assert isinstance(
                df.seqs.seq.back_translate(codon_map=codon_table_to_codon_map(human_codon_table)),
                pd.Series,
            )


# ---- test list of seqs of various types ---------
seqs_list_to_numpy = [
    (seqs, "nt", 7),
    (seqs_aa, "aa", 23),
    pytest.param(
        seqs_mixed,
        None,
        None,
        marks=pytest.mark.xfail(reason="Not a homogeneous list of SeqLikes."),
    ),
]


@pytest.mark.parametrize("seqs, _type, n_bases", seqs_list_to_numpy)
def test_to_onehot(seqs, _type, n_bases):
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": seqs})

    assert isinstance(df.seqs.seq.to_onehot(), np.ndarray)

    num_seqs, len_seq, num_bases = df.seqs.seq.to_onehot().shape
    assert num_seqs == len(df)
    assert len_seq == len(df.seqs.iloc[0])
    assert num_bases == n_bases
    assert df.seqs.seq.to_onehot().max() == 1
    assert df.seqs.seq.to_onehot().min() == 0


@pytest.mark.parametrize("seqs, _type, n_bases", seqs_list_to_numpy)
def test_to_index(seqs, _type, n_bases):
    # TODO: Docstring needed for test intent.
    df = pd.DataFrame({"seqs": seqs})

    assert isinstance(df.seqs.seq.to_index(), np.ndarray)

    num_seqs, len_seq = df.seqs.seq.to_index().shape
    assert num_seqs == len(df)
    assert len_seq == len(df.seqs.iloc[0])
    assert df.seqs.seq.to_index().max() <= n_bases
    assert df.seqs.seq.to_index().min() >= 0
