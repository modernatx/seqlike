import pytest
from seqlike import SeqLike

seq = "ATGCATGCATGCATGCAT"


def test_len():
    assert len(SeqLike(seq, "dna").nt()) == 18
    assert len(SeqLike(seq, "dna").aa()) == 6
    assert len(SeqLike(seq, "dna")) == 18


def test_contain():
    assert ("AUG" in SeqLike(seq, "dna").nt()) == False
    assert ("CAT" in SeqLike(seq, "dna").nt()) == True
    assert ("CAT" in SeqLike(seq, "dna").aa()) == False
    assert ("HAC" in SeqLike(seq, "dna").aa()) == True
    assert ("AUG" in SeqLike(seq, "dna")) == False


def test_iter():
    iter_list = []
    for index, nucleotide in enumerate(SeqLike(seq, "dna").nt()):
        if 5 < index < 10:
            iter_list.append(nucleotide)

    assert iter_list == ["G", "C", "A", "T"]

    for index, nucleotide in enumerate(SeqLike(seq, "dna")):
        print(nucleotide)


def test_slicing():
    s = SeqLike(seq, "dna")[5]
    assert isinstance(s, SeqLike)
    assert s._type == "NT"
    assert str(s) == str(SeqLike("T", seq_type="dna"))

    s = SeqLike(seq, "dna").nt()[5]
    assert isinstance(s, SeqLike)
    assert s._type == "NT"
    assert str(s) == str(SeqLike("T", seq_type="dna"))

    s = SeqLike(seq, "dna").aa()[3]
    assert isinstance(s, SeqLike)
    assert s._type == "AA"
    assert str(s) == str(SeqLike("C", seq_type="aa"))

    s = SeqLike(seq, "dna").nt()[1:5]
    # if we slice to a non-multiple of 3, we shouldn't be able to
    # convert to aa and if we try, it should raise a type error
    with pytest.raises(TypeError):
        s.aa()
    assert isinstance(s, SeqLike)
    assert s._type == "NT"
    assert str(s) == str(SeqLike("TGCA", seq_type="dna"))
    assert s._aa_record == None

    # if we slice a multiple of 3, we should be able to translate
    s = SeqLike(seq, "dna")[:3].aa()
    assert isinstance(s, SeqLike)
    assert s._type == "AA"
    assert str(s) == str(SeqLike("M", seq_type="aa"))
    assert s._type == "AA"

    # if we slice an aa, then we should get the NT (if it was there)
    s = SeqLike(seq, "dna").aa()[1:5]
    assert isinstance(s, SeqLike)
    assert s._type == "AA"
    assert str(s) == str(SeqLike("HACM", seq_type="aa"))
    assert len(s._nt_record) == 12

    # but if we never had a NT associated, then we shouldn't expect
    # one in the sliced result
    s = SeqLike("FEFEFEE", "aa")[:3]
    assert s._type == "AA"
    assert s._nt_record == None
