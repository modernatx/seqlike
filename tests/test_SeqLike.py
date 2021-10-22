import os
import sys
from copy import deepcopy
import tempfile

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import pytest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from seqlike import SeqLike

from seqlike.codon_tables import yeast_codon_table, ecoli_codon_table, codon_table_to_codon_map
from seqlike.SeqLike import NT, AA, STANDARD_NT, STANDARD_AA

from . import test_path


from hypothesis import given, assume
from hypothesis.strategies import composite, text, integers, sampled_from


ALPHABETS = [NT, AA, STANDARD_NT, STANDARD_AA]
NT_TYPES = ["NT", "DNA", "RNA", "dna", "nt", "rna"]
AA_TYPES = ["AA", "aa"]
SEQ_TYPES = NT_TYPES + AA_TYPES


@composite
def string_sequences(draw, alphabet: str = None, min_size: int = 1) -> tuple:
    """Composite Hypothesis strategy to generate strings from one of the alphabets.

    This test handles the cases where seq_type, sequence, and alphabet are all set.
    """
    max_size = draw(integers(min_value=min_size, max_value=1000))
    if not alphabet:
        alphabet = draw(sampled_from(ALPHABETS))
    sequence = draw(text(alphabet=alphabet, min_size=min_size, max_size=max_size))
    if alphabet in [NT, STANDARD_NT]:
        seq_type = draw(sampled_from(NT_TYPES))
    elif alphabet in [AA, STANDARD_AA]:
        seq_type = draw(sampled_from(AA_TYPES))
    return sequence, seq_type, alphabet


@given(string_sequences())
def test_init_hyp(sequence_type_and_alphabet: tuple):
    """Test initialization of SeqLike objects from strings."""
    sequence, seq_type, alphabet = sequence_type_and_alphabet
    seq = SeqLike(sequence, alphabet=alphabet, seq_type=seq_type)
    assert seq.alphabet == alphabet

    if alphabet == AA or alphabet == STANDARD_AA:
        assert seq._type == "AA"
        assert seq._nt_record is None
    elif alphabet == NT or alphabet == STANDARD_NT:
        assert seq._type == "NT"
        assert seq._aa_record is None


# 26 July 2021
# We need a test that breaks the constructor.


def test_init_fastq():
    seq0 = SeqLike(SeqIO.read(test_path / f"sanger.ab1", "abi"), "nt")
    assert list(seq0.letter_annotations.keys()) == ["phred_quality", "seqnums"]
    with tempfile.NamedTemporaryFile(mode="w") as tempf:
        SeqIO.write(seq0.to_seqrecord(), tempf, "fastq")


def test_init_custom_alphabet():
    """Test that custom alphabets are passed right through."""
    sequence = "ATGC"
    seq = SeqLike(sequence, "nt", alphabet="*-ACDEFG")
    assert seq.alphabet == "*-ACDEFG"


@given(string_sequences())
def test_init_from_seqrecord(sequence_type_and_alphabet):
    """Test that initialization from SeqRecord works."""
    sequence, seq_type, alphabet = sequence_type_and_alphabet
    details = {"id": "i1", "name": "n1", "description": "desc1"}
    seqrec = SeqLike(sequence, seq_type=seq_type, alphabet=alphabet).to_seqrecord(**details)
    seq = SeqLike(seqrec, seq_type=seq_type)
    assert seq.id == "i1"
    assert seq.name == "n1"
    assert seq.description == "desc1"


@given(string_sequences())
def test_init_from_seqrecord_copied(sequence_type_and_alphabet):
    """Test that initialization from SeqRecords return deep copies."""
    sequence, seq_type, alphabet = sequence_type_and_alphabet
    seqrec = SeqLike(sequence, alphabet=alphabet, seq_type=seq_type, id="id0").to_seqrecord()

    seq0 = SeqLike(seqrec, seq_type=seq_type, id=f"1:{seqrec.id}")
    assert seqrec.id == "id0"
    assert seq0.id == "1:id0"
    seq1 = SeqLike(seq0, seq_type=seq_type, id=f"2:{seq0.id}")
    assert seq0.id == "1:id0"
    assert seq1.id == "2:1:id0"


def test_SeqLike_interconversion():
    """Test that SeqLike interconversion works properly.

    Suggested TODO: Split out SeqLike interconversion into a parametrized interconversion
    based on target types.
    """
    seq = "TCGCACACTGCA"

    a1 = SeqLike(seq, "dna").nt().to_str()
    a2 = SeqLike(seq, "dna").aa().to_str()
    assert isinstance(a1, str)
    assert isinstance(a2, str)

    b1 = SeqLike(seq, "dna").nt().to_seq()
    b2 = SeqLike(seq, "dna").aa().to_seq()
    assert isinstance(b1, Seq)
    assert isinstance(b2, Seq)

    c1 = SeqLike(seq, "dna").nt().to_seqrecord()
    c2 = SeqLike(seq, "dna").aa().to_seqrecord()
    assert isinstance(c1, SeqRecord)
    assert isinstance(c2, SeqRecord)

    seq1 = SeqLike(seq, "dna").nt()
    seq2 = SeqLike(seq, "dna").aa()
    d1 = seq1.to_onehot()
    d2 = seq2.to_onehot()
    assert isinstance(d1, np.ndarray)
    assert isinstance(d2, np.ndarray)

    e1 = SeqLike(d1, "dna", alphabet=seq1.alphabet).to_str()
    e2 = SeqLike(d2, "aa", alphabet=seq2.alphabet).to_str()
    assert e1 == a1
    assert e2 == a2
    # when interconverting, should the letter_annotations be empty?
    assert seq2.letter_annotations == {}
    seq3 = SeqLike(seq2, "dna")
    seqnums = [str(i + 1) for i in range(len(seq3))]
    assert seq3.letter_annotations["seqnums"] == seqnums
    assert seq3[:2].letter_annotations["seqnums"] == seqnums[:2]


CODON_TABLES = {"ecoli_k12": yeast_codon_table, "Kazusa_yeast": yeast_codon_table}


@composite
def pick_codon_map(draw):
    # TODO: Add docstrings here.
    table = draw(sampled_from(list(CODON_TABLES.keys())))
    return CODON_TABLES[table]


@given(string_sequences(alphabet=STANDARD_AA, min_size=10))
def test_backtranslation_without_codon_map(sequence_type_and_alphabet):
    """Test that back-translation without codon map set raises AttributeError."""
    sequence, seq_type, _ = sequence_type_and_alphabet
    s = SeqLike(
        sequence,
        id="test1",
        seq_type=seq_type,
        description="desc1",
        annotations={"a1": "test1"},
    )
    with pytest.raises(AttributeError):
        s.back_translate()


@given(pick_codon_map(), string_sequences(alphabet=STANDARD_AA, min_size=10))
def test_backtranslation_overriding_codon_map(codon_table, sequence_type_and_alphabet):
    """Test that back-translation works when we override the codon map."""
    sequence, seq_type, _ = sequence_type_and_alphabet
    codon_map = codon_table_to_codon_map(codon_table, deterministic=False)
    s = SeqLike(
        sequence,
        id="test1",
        seq_type=seq_type,
        description="desc1",
        annotations={"a1": "test1"},
    )
    s.codon_map = codon_map
    assert_back_translate_properties(s)


def assert_back_translate_properties(s):
    """Test for back_translate properties.

    This collection of assertions repeatedly comes up
    so we have refactored them out into a separate function.
    """
    assert isinstance(s.back_translate(), SeqLike)
    # retains SeqLike metadata like `id` and `description`
    assert s.back_translate().id == s.id
    assert s.back_translate().description == s.description
    assert s.back_translate().annotations == s.annotations
    assert s.back_translate().dbxrefs == s.dbxrefs
    # but letter-specific metadata is changed
    assert s.back_translate().letter_annotations != s.letter_annotations

    # works, overrides self.codon_map
    # assert isinstance(s2.back_translate(codon_map=codon_map2), SeqLike)

    seqs = [s.back_translate().nt().to_str() for _ in range(10)]
    assert len(seqs) == 10
    assert len(set(seqs)) <= len(seqs)

    # DNA sequences are also valid AA sequences,
    # and we can only backtranslate a AA sequence.
    with pytest.raises(TypeError):
        s.nt().back_translate()


@given(pick_codon_map(), string_sequences(alphabet=STANDARD_AA, min_size=10))
def test_backtranslation(codon_table, sequence_type_and_alphabet):
    """Test that back-translation works."""
    sequence, seq_type, _ = sequence_type_and_alphabet
    codon_map = codon_table_to_codon_map(codon_table, deterministic=False)
    s = SeqLike(
        sequence,
        id="test1",
        seq_type=seq_type,
        description="desc1",
        annotations={"a1": "test1"},
        codon_map=codon_map,
    )


def test_upper():
    # TODO: Add docstring
    # TODO: Parametrize with hypothesis.
    s1 = SeqLike("TCgCAcActgcA", seq_type="nt")
    s2 = SeqLike("GEgdatygKLTlkfiCTT", seq_type="aa")
    assert isinstance(s1.upper(), SeqLike)
    assert isinstance(s2.upper(), SeqLike)
    assert str(s1.upper()) == str(s1).upper()
    assert str(s2.upper()) == str(s2).upper()


def test_count():
    # TODO: Add docstring
    # TODO: Parametrize with hypothesis.
    s1 = SeqLike("TCGCACACTGCA", seq_type="nt")
    s2 = SeqLike("GEGDRTYNKLTLIFQCTT", seq_type="aa")
    assert isinstance(s1.count("T"), int)
    assert isinstance(s2.count("Y"), int)
    for letter in s1:
        assert s1.count(letter) == str(s1).count(letter)
    for letter in s2:
        assert s2.count(letter) == str(s2).count(letter)


def test_find():
    # TODO: Add docstring
    # TODO: Parametrize with hypothesis.
    s1 = SeqLike("TCGCACACTGCATCGCACACTGCA", seq_type="dna")

    sub1 = SeqLike("ACTG", seq_type="dna")
    found = s1.find(sub1)
    print(found)
    assert isinstance(found, int)
    assert found == 6
    assert s1.find(sub1) == 6
    assert s1.find(sub1, 1) == 6
    assert s1.find(sub1, 7) == 18
    assert s1.find("CCCCCCC") == -1

    for substr in ["ACTG", "TCGCACA"]:
        assert s1.find(SeqLike(substr, "dna")) == s1.find(substr)
    for start, stop in [(0, 5), (8, 11), (5, 15)]:
        assert s1.find(s1[start:stop]) == s1.find(str(s1)[start:stop])


def test_reverse_complement():
    # TODO: Add docstring
    # TODO: Parametrize with hypothesis.
    s1 = SeqLike("TTTCACACTGCA", seq_type="nt")
    assert isinstance(s1.reverse_complement(), SeqLike)
    assert str(s1.reverse_complement()) == str(s1.to_seq().reverse_complement())
    assert s1.reverse_complement().id == s1.id
    assert s1.reverse_complement().name == s1.name
    assert s1.reverse_complement().annotations["reversed"]

    assert str(s1.reverse_complement().reverse_complement()) == str(s1)
    assert not s1.reverse_complement().reverse_complement().annotations["reversed"]
    with pytest.raises(ValueError):
        s2 = SeqLike("GEGDATYGKLTLKFICTT", seq_type="aa")
        s2.reverse_complement()

    ## .aa() is of the reverse complemented sequence
    s3 = SeqLike("TGCAGTGTGAAA", seq_type="nt")
    assert str(s3.reverse_complement()) == str(s1)
    assert str(s3.reverse_complement().aa()) == str(s1.aa())


def test_apply():
    # TODO: Add docstring
    # TODO: Parametrize with hypothesis.
    codon_map = codon_table_to_codon_map(ecoli_codon_table, deterministic=True)

    s = SeqLike("GEGDATYGKLTLKFICTT", seq_type="aa", codon_map=codon_map)
    s2 = SeqLike("GEGDATYGKLTLKFICTT", seq_type="aa")

    # test an apply that returns a non-SeqLike
    assert s.apply(count_residue) == 4

    # test an apply that returns a SeqLike
    s_split = s.apply(split_seq, split_location=5)
    assert len(s_split) == 5
    assert isinstance(s_split, SeqLike)

    assert s.apply(split_seq, split_location=5).back_translate().codon_map is not None
    assert s_split.codon_map is not None


def test_ungap():
    seqstr = "GEGDATYGK--LTLKFICTT"
    s = SeqLike(seqstr, seq_type="aa")
    assert s.ungap().to_str() == seqstr.replace("-", "")


def test_seq_num_to_idx():
    s = SeqLike("N" * 20, seq_type="nt")
    assert s.seq_num_to_idx(["2", "3", "4", "5"]) == [slice(1, 5)]
    assert s.seq_num_to_idx(["2", "3", "4", "5", "1", "2"]) == [
        slice(1, 5),
        slice(0, 2),
    ]
    assert s.seq_num_to_idx(["2", "3", "4", "5", "1", "3"]) == [slice(1, 5), 0, 2]


@pytest.mark.parametrize("seqstr, seq_type", [("TCGCACACTGCA", "nt"), ("GEGDATYGKLTLKFICTT", "aa")])
def test_slice(seqstr, seq_type):
    s = SeqLike(seqstr, seq_type=seq_type)
    for seqnums in ["2", "3", "4", "5"], ["1", "10", "12"], ["10", "11", "12", "1"]:
        assert s.slice(seqnums).to_str() == "".join(seqstr[int(n) - 1] for n in seqnums)
        assert s.slice(seqnums).id == s.id


def test_slice_with_insertion_codes():
    # test numbering with insertion codes, e.g., antibody numbering
    seqstr = "GEGDATY"
    seqnums = "99,100,100A,100B,100C,101,102".split(",")
    s = SeqLike(seqstr, seq_type="aa")
    s.letter_annotations["seqnums"] = seqnums
    assert s.slice("99,100,101,102".split(",")).to_str() == "GETY"
    assert s.slice("100,100A,100B".split(",")).to_str() == "EGD"
    assert s.slice("100C,99,101".split(",")).to_str() == "AGT"


def test__getattr__():
    for seqstr, seq_type in [("TCGCACACTGCA", "nt"), ("GEGDATYGKLTLKFICTT", "aa")]:
        seqid = "test1"
        name = "test2"
        description = "test3"
        annotations = {"reversed": True}
        letter_annotations = {"seqnums": [str(i + 1) for i in range(len(seqstr))]}
        seqrec = SeqRecord(
            Seq(seqstr),
            id=seqid,
            name=name,
            description=description,
            annotations=annotations,
            letter_annotations=letter_annotations,
        )
        seq = SeqLike(seqrec, seq_type)
        assert str(seq) == seqstr
        assert seq.id == seqrec.id == seqid
        assert seq.name == seqrec.name == name
        assert seq.description == seqrec.description == description
        assert seq.annotations == seqrec.annotations == annotations
        assert seq.letter_annotations == seqrec.letter_annotations == letter_annotations


def test__setattr__():
    for seqstr, seq_type in [("TCGCACACTGCA", "nt"), ("GEGDATYGKLTLKFICTT", "aa")]:
        seqid = "test1"
        name = "test2"
        description = "test3"
        annotations = {"reversed": True}
        # this is intentionally zero-index, as SeqLike automatically generates 1-index seqnum annotations
        letter_annotations = {"seqnums": [str(i) for i in range(len(seqstr))]}
        seqrec = SeqRecord(
            Seq(seqstr),
            id=seqid,
            name=name,
            description=description,
            annotations=annotations,
            letter_annotations=letter_annotations,
        )
        # set values during initialization
        seq = SeqLike(
            seqstr,
            seq_type=seq_type,
            id=seqid,
            name=name,
            description=description,
            annotations=annotations,
            letter_annotations=letter_annotations,
        )
        assert str(seq) == seqstr
        assert seq.id == seq._seqrecord.id == seqid
        assert seq.name == seq._seqrecord.name == name
        assert seq.description == seq._seqrecord.description == description
        assert seq.annotations == seq._seqrecord.annotations == annotations
        assert seq.letter_annotations == seq._seqrecord.letter_annotations == letter_annotations

        # set values after init
        seqid = "new test1"
        name = "new test2"
        description = "new test3"
        annotations = {"reversed": True, "annot2": 2}
        letter_annotations = {
            "seqnums": [str(i) for i in range(len(seqstr))],
            "sec_str": ["S"] * len(seqstr),
        }
        seq.id = seqid
        seq.name = name
        seq.description = description
        seq.annotations = annotations
        seq.letter_annotations = letter_annotations
        assert seq.id == seq._seqrecord.id == seqid
        assert seq.name == seq._seqrecord.name == name
        assert seq.description == seq._seqrecord.description == description
        assert seq.annotations == seq._seqrecord.annotations == annotations
        assert seq.letter_annotations == seq._seqrecord.letter_annotations == letter_annotations
        # should be different from original values
        assert seq.id != seqrec.id
        assert seq.name != seqrec.name
        assert seq.description != seqrec.description


def test__dir__():
    for seqstr, seq_type in [("TCGCACACTGCA", "nt"), ("GEGDATYGKLTLKFICTT", "aa")]:
        seqid = "test1"
        name = "test2"
        description = "test3"
        annotations = {"reversed": True}
        seq = SeqLike(
            seqstr,
            seq_type=seq_type,
            id=seqid,
            name=name,
            description=description,
            annotations=annotations,
        )
        assert isinstance(dir(seq), list)
        assert len(dir(seq)) == 68
        assert set(dir(SeqLike)).issubset(dir(seq))
        assert set(["annotations", "description", "name", "id", "letter_annotations"]).issubset(dir(seq))


def test__getitem__():
    # TODO: Parametrize this test.
    for seqstr, seq_type in [("TCGCACACTGCA", "nt"), ("GEGDATYGKLTLKFICTT", "aa")]:
        seq = SeqLike(seqstr, seq_type)
        assert type(seq) == SeqLike
        assert str(seq[0]) == seqstr[0]  # "T"
        assert str(seq[:1]) == seqstr[:1]  # "T"
        assert str(seq[:3]) == seqstr[:3]  # "TCG"
        assert str(seq[2:6]) == seqstr[2:6]  # "GCAC"
        assert str(seq[-5:]) == seqstr[-5:]  # "CTGCA"
        assert str(seq[-5:-1]) == seqstr[-5:-1]  # "CTGC"
        assert str(seq[:0]) == seqstr[:0]  # ""
        assert str(seq[::-1]) == seqstr[::-1]  # "ACGTCACACGCT"

        seqnums = seq.letter_annotations["seqnums"]
        assert seq[0].letter_annotations["seqnums"] == [seqnums[0]]
        assert seq[:1].letter_annotations["seqnums"] == seqnums[:1]
        assert seq[:3].letter_annotations["seqnums"] == seqnums[:3]
        assert seq[2:6].letter_annotations["seqnums"] == seqnums[2:6]
        assert seq[-5:].letter_annotations["seqnums"] == seqnums[-5:]
        assert seq[-5:-1].letter_annotations["seqnums"] == seqnums[-5:-1]
        assert seq[:0].letter_annotations["seqnums"] == seqnums[:0]
        assert seq[::-1].letter_annotations["seqnums"] == seqnums[::-1]


def test__add__():
    # TODO: Parametrize this test.
    """Test for __add__ method."""
    for seqstr, seq_type in [("TCGCACACTGCA", "nt"), ("GEGDATYGKLTLKFICTT", "aa")]:
        annotations = {"a1": "test1"}
        s1 = SeqLike(seqstr[:5], seq_type=seq_type, description="test", annotations=annotations)
        s2 = SeqLike(seqstr[5:], seq_type=seq_type, annotations=annotations)
        s3 = SeqLike(seqstr[5:], seq_type=seq_type)
        assert (s1 + s2).to_str() == seqstr
        assert (s1 + s2).id == s1.id
        assert (s1 + s2).description == s1.description == "test"
        # when annotations between two SeqLike match, the concatenated object has same
        assert (s1 + s2).annotations == s1.annotations == annotations
        # if not, the concatenated object does not (to match behavior of SeqRecord)
        assert (s1 + s3).annotations != s1.annotations
        assert (s1 + s2).dbxrefs == s1.dbxrefs
        # for __radd__
        assert sum([s1, s2]).to_str() == seqstr
        assert sum([s1, s2]).id == s1.id
        assert sum([s1, s2]).description == s1.description == "test"
        assert sum([s1, s2, s1]).id == s1.id
        assert sum([s1]).to_str() == str(s1)
        # for __iadd__
        s12 = SeqLike(s1, seq_type=seq_type)
        s12 += s2
        # adding a string
        s12s = SeqLike(s1, seq_type=seq_type)
        s12s += str(s2)
        assert s12.to_str() == s12s.to_str() == seqstr
        assert s12.id == s12s.id == s1.id
        assert s12.description == s12s.description == s1.description == "test"


def test__deepcopy__():
    for seqstr, seq_type in [("TCGCACACTGCA", "nt"), ("GEGDATYGKLTLKFICTT", "aa")]:
        s = SeqLike(seqstr, seq_type=seq_type, description="test")
        assert deepcopy(s).to_str() == seqstr
        assert deepcopy(s).id == s.id
        assert deepcopy(s).description == s.description == "test"
        # deepcopy returns a new object, so it should not match the original
        assert deepcopy(s) != s


def split_seq(seq, split_location):
    return seq[:split_location]


def count_residue(seq, residue="T"):
    return seq.to_str().count(residue)
