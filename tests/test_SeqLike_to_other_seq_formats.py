from types import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from seqlike import SeqLike
from seqlike.alphabets import is_NT, is_AA
import pytest


# A list of sequences to test for .nt() call-
# DNA, RNA, Protein and non-biological. Since .nt() call for a protein input requires back translation, pytest for it will fail.
sequence_list_nt = [
    ("ATGCGCGCG", "ATGCGCGCG", True, "dna"),
    ("AUGGUGCUGUUU", "AUGGUGCUGUUU", True, "rna"),
    pytest.param(
        "LAIGHATY",
        "LAIGHATY",
        True,
        "aa",
        marks=pytest.mark.xfail(reason=".nt() call expects to back translate AA sequence."),
    ),
    pytest.param(
        "#^@*)!&$&$)_#_",
        "#^@*)!&$&$)_#_",
        True,
        None,
        marks=pytest.mark.xfail(reason="Not a biological sequence."),
    ),
]


# A list of sequences to test for .aa() call- DNA, RNA, Protein and non-biological.
sequence_list_aa = [
    ("ATGCGCGCG", "MRA", True, "dna"),
    ("AUGGUGCUGUUU", "MVLF", True, "rna"),
    ("LAIGHATY", "LAIGHATY", True, "aa"),
    pytest.param(
        "#^@*)!&$&$)_#_",
        "#^@*)!&$&$)_#_",
        True,
        None,
        marks=pytest.mark.xfail(reason="Not a biological sequence."),
    ),
]


# String input to string
@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_nt)
def test_str_nt_str(sequence, expected_sequence, expected_boolean, seq_type):
    nt_seq = SeqLike(sequence, seq_type=seq_type).nt().to_str()

    assert nt_seq == expected_sequence
    assert isinstance(nt_seq, str), "string conversion failed, %s" % type(nt_seq)
    assert is_NT(nt_seq) == expected_boolean


@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_aa)
def test_str_aa_str(sequence, expected_sequence, expected_boolean, seq_type):
    aa_seq = SeqLike(sequence, seq_type=seq_type).aa().to_str()

    assert aa_seq == expected_sequence
    assert isinstance(aa_seq, str), "string conversion failed, %s" % type(aa_seq)
    assert is_AA(aa_seq) == expected_boolean


# String input to Seq object
@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_nt)
def test_str_nt_seq(sequence, expected_sequence, expected_boolean, seq_type):
    nt_seq = SeqLike(sequence, seq_type=seq_type).nt().to_seq()
    nt_seq_seq = str(nt_seq)

    assert nt_seq_seq == expected_sequence
    assert isinstance(nt_seq, Seq), "string conversion failed, %s" % type(nt_seq)
    assert is_NT(nt_seq_seq) == expected_boolean


@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_aa)
def test_str_aa_seq(sequence, expected_sequence, expected_boolean, seq_type):
    aa_seq = SeqLike(sequence, seq_type=seq_type).aa().to_seq()
    aa_seq_seq = str(aa_seq)

    assert aa_seq_seq == expected_sequence
    assert isinstance(aa_seq, Seq), "string conversion failed, %s" % type(aa_seq)
    assert is_AA(aa_seq_seq) == expected_boolean


# String input to SeqRecord object
@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_nt)
def test_str_nt_seqrecord(sequence, expected_sequence, expected_boolean, seq_type):
    nt_seq = SeqLike(sequence, seq_type=seq_type).nt().to_seqrecord()
    nt_seq_seq = str(nt_seq.seq)

    assert nt_seq_seq == expected_sequence
    assert isinstance(nt_seq, SeqRecord), "string conversion failed, %s" % type(nt_seq)
    assert is_NT(nt_seq_seq) == expected_boolean


@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_aa)
def test_str_aa_seqrecord(sequence, expected_sequence, expected_boolean, seq_type):
    aa_seq = SeqLike(sequence, seq_type=seq_type).aa().to_seqrecord()
    aa_seq_seq = str(aa_seq.seq)

    assert aa_seq_seq == expected_sequence
    assert isinstance(aa_seq, SeqRecord), "string conversion failed, %s" % type(aa_seq)
    assert is_AA(aa_seq_seq) == expected_boolean


# A list of Seq objects to test for .nt() call-
# DNA, RNA, Protein and non-biological. Since .nt() call for a protein input requires back translation, pytest for it will fail.
sequence_list_nt = [
    (Seq("ATGCGCGCG"), "ATGCGCGCG", True, "dna"),
    (Seq("AUGGUGCUGUUU"), "AUGGUGCUGUUU", True, "rna"),
    pytest.param(
        Seq("LAIGHATY"),
        "LAIGHATY",
        True,
        "aa",
        marks=pytest.mark.xfail(reason=".nt() call expects to back translate input AA sequence."),
    ),
    pytest.param(
        Seq("#^@*)!&$&$)_#_"),
        "#^@*)!&$&$)_#_",
        True,
        None,
        marks=pytest.mark.xfail(reason="Not a biological sequence."),
    ),
]


# A list of Seq objects to test for .aa() call- DNA, RNA, Protein and non-biological.
sequence_list_aa = [
    (Seq("ATGCGCGCG"), "MRA", True, "dna"),
    (Seq("AUGGUGCUGUUU"), "MVLF", True, "rna"),
    (Seq("LAIGHATY"), "LAIGHATY", True, "aa"),
    (
        Seq(
            "LAIGHATY",
        ),
        "LAIGHATY",
        True,
        "aa",
    ),
    pytest.param(
        Seq("#^@*)!&$&$)_#_"),
        "#^@*)!&$&$)_#_",
        True,
        None,
        marks=pytest.mark.xfail(reason="Not a biological sequence."),
    ),
]


# Seq object input to string
@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_nt)
def test_seq_nt_str(sequence, expected_sequence, expected_boolean, seq_type):
    nt_seq = SeqLike(sequence, seq_type=seq_type).nt().to_str()

    assert nt_seq == expected_sequence
    assert isinstance(nt_seq, str), "string conversion failed, %s" % type(nt_seq)
    assert is_NT(nt_seq) == expected_boolean


@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_aa)
def test_seq_aa_str(sequence, expected_sequence, expected_boolean, seq_type):
    aa_seq = SeqLike(sequence, seq_type=seq_type).aa().to_str()

    assert aa_seq == expected_sequence
    assert isinstance(aa_seq, str), "string conversion failed, %s" % type(aa_seq)
    assert is_AA(aa_seq) == expected_boolean


# Seq object input to Seq object
@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_nt)
def test_seq_nt_seq(sequence, expected_sequence, expected_boolean, seq_type):
    nt_seq = SeqLike(sequence, seq_type=seq_type).nt().to_seq()
    nt_seq_seq = str(nt_seq)

    assert nt_seq_seq == expected_sequence
    assert isinstance(nt_seq, Seq), "string conversion failed, %s" % type(nt_seq)
    assert is_NT(nt_seq_seq) == expected_boolean


@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_aa)
def test_seq_aa_seq(sequence, expected_sequence, expected_boolean, seq_type):
    aa_seq = SeqLike(sequence, seq_type=seq_type).aa().to_seq()
    aa_seq_seq = str(aa_seq)

    assert aa_seq_seq == expected_sequence
    assert isinstance(aa_seq, Seq), "string conversion failed, %s" % type(aa_seq)
    assert is_AA(aa_seq_seq) == expected_boolean


# Seq object input to SeqRecord object
@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_nt)
def test_seq_nt_seqrecord(sequence, expected_sequence, expected_boolean, seq_type):
    nt_seq = SeqLike(sequence, seq_type=seq_type).nt().to_seqrecord()
    nt_seq_seq = str(nt_seq.seq)

    assert nt_seq_seq == expected_sequence
    assert isinstance(nt_seq, SeqRecord), "string conversion failed, %s" % type(nt_seq)
    assert is_NT(nt_seq_seq) == expected_boolean


@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_aa)
def test_seq_aa_seqrecord(sequence, expected_sequence, expected_boolean, seq_type):
    aa_seq = SeqLike(sequence, seq_type=seq_type).aa().to_seqrecord()
    aa_seq_seq = str(aa_seq.seq)

    assert aa_seq_seq == expected_sequence
    assert isinstance(aa_seq, SeqRecord), "string conversion failed, %s" % type(aa_seq)
    assert is_AA(aa_seq_seq) == expected_boolean


# A list of Seq objects to test for .nt() call-
# DNA, RNA, Protein and non-biological. Since .nt() call for a protein input requires back translation, pytest for it will fail.
sequence_list_nt = [
    (SeqRecord(Seq("ATGCGCGCG")), "ATGCGCGCG", True, "dna"),
    (SeqRecord(Seq("AUGGUGCUGUUU")), "AUGGUGCUGUUU", True, "rna"),
    pytest.param(
        SeqRecord(Seq("LAIGHATY")),
        "LAIGHATY",
        True,
        "aa",
        marks=pytest.mark.xfail(reason=".nt() call expects to back translate input AA sequence."),
    ),
    pytest.param(
        SeqRecord(Seq("#^@*)!&$&$)_#_")),
        "#^@*)!&$&$)_#_",
        True,
        None,
        marks=pytest.mark.xfail(reason="Not a biological sequence."),
    ),
]


# A list of Seq objects to test for .aa() call- DNA, RNA, Protein and non-biological.
sequence_list_aa = [
    (SeqRecord(Seq("ATGCGCGCG")), "MRA", True, "dna"),
    (SeqRecord(Seq("AUGGUGCUGUUU")), "MVLF", True, "rna"),
    (SeqRecord(Seq("LAIGHATY")), "LAIGHATY", True, "aa"),
    pytest.param(
        SeqRecord(Seq("#^@*)!&$&$)_#_")),
        "#^@*)!&$&$)_#_",
        True,
        None,
        marks=pytest.mark.xfail(reason="Not a biological sequence."),
    ),
]


# SeqRecord input to string
@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_nt)
def test_seqrecord_nt_str(sequence, expected_sequence, expected_boolean, seq_type):
    nt_seq = SeqLike(sequence, seq_type=seq_type).nt().to_str()

    assert nt_seq == expected_sequence
    assert isinstance(nt_seq, str), "string conversion failed, %s" % type(nt_seq)
    assert is_NT(nt_seq) == expected_boolean


@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_aa)
def test_seqrecord_aa_str(sequence, expected_sequence, expected_boolean, seq_type):
    aa_seq = SeqLike(sequence, seq_type=seq_type).aa().to_str()

    assert aa_seq == expected_sequence
    assert isinstance(aa_seq, str), "string conversion failed, %s" % type(aa_seq)
    assert is_AA(aa_seq) == expected_boolean


# Seqrecord input to Seq object
@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_nt)
def test_seqrecord_nt_seq(sequence, expected_sequence, expected_boolean, seq_type):
    nt_seq = SeqLike(sequence, seq_type=seq_type).nt().to_seq()
    nt_seq_seq = str(nt_seq)

    assert nt_seq_seq == expected_sequence
    assert isinstance(nt_seq, Seq), "string conversion failed, %s" % type(nt_seq)
    assert is_NT(nt_seq_seq) == expected_boolean


@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_aa)
def test_seqrecord_aa_seq(sequence, expected_sequence, expected_boolean, seq_type):
    aa_seq = SeqLike(sequence, seq_type=seq_type).aa().to_seq()
    aa_seq_seq = str(aa_seq)

    assert aa_seq_seq == expected_sequence
    assert isinstance(aa_seq, Seq), "string conversion failed, %s" % type(aa_seq)
    assert is_AA(aa_seq_seq) == expected_boolean


# SeqRecord input to SeqRecord object
@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_nt)
def test_seqrecord_nt_seqrecord(sequence, expected_sequence, expected_boolean, seq_type):
    nt_seq = SeqLike(sequence, seq_type=seq_type).nt().to_seqrecord()
    nt_seq_seq = str(nt_seq.seq)

    assert nt_seq_seq == expected_sequence
    assert isinstance(nt_seq, SeqRecord), "string conversion failed, %s" % type(nt_seq)
    assert is_NT(nt_seq_seq) == expected_boolean


@pytest.mark.parametrize("sequence, expected_sequence, expected_boolean, seq_type", sequence_list_aa)
def test_seqrecord_aa_seqrecord(sequence, expected_sequence, expected_boolean, seq_type):
    aa_seq = SeqLike(sequence, seq_type=seq_type).aa().to_seqrecord()
    aa_seq_seq = str(aa_seq.seq)

    assert aa_seq_seq == expected_sequence
    assert isinstance(aa_seq, SeqRecord), "string conversion failed, %s" % type(aa_seq)
    assert is_AA(aa_seq_seq) == expected_boolean
