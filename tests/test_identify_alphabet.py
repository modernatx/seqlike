import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from seqlike.SeqLike import SeqLike
from seqlike.alphabets import AA, is_AA, is_NT


def test_rna():
    seqstr = "AUGC"  # can only be RNA or AA (with extented alphabet)
    for seq in [
        seqstr,
        seqstr.lower(),
        Seq(seqstr),
        SeqRecord(Seq(seqstr)),
        SeqLike(seqstr, seq_type="nt"),
    ]:
        assert is_NT(seq)
        assert is_AA(seq)


def test_dna():
    seqstr = "ATCG"
    for seq in [seqstr, seqstr.lower(), Seq(seqstr), SeqRecord(Seq(seqstr)), SeqLike(seqstr, seq_type="nt")]:
        assert is_NT(seq)
        assert is_AA(seq)


def test_protein():
    seqstr = "GALMFACT"
    for seq in [seqstr, seqstr.lower(), Seq(seqstr), SeqRecord(Seq(seqstr)), SeqLike(seqstr, seq_type="aa")]:
        assert is_AA(seq)
        assert not is_NT(seq)


def test_non_biological():
    seq = "$$$"
    assert not is_NT(seq)
    assert not is_AA(seq)
