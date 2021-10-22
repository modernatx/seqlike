import os
import sys
from copy import deepcopy
import tempfile

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import pytest

from seqlike.SeqLike import NT, AA, STANDARD_NT, STANDARD_AA, SeqLike, aaSeqLike, ntSeqLike, is_AA, is_NT

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
def test_NTSeqLike(sequence_type_and_alphabet):
    """Test initialization of SeqLike objects from strings."""
    sequence, seq_type, alphabet = sequence_type_and_alphabet
    seq = SeqLike(sequence, alphabet=alphabet, seq_type=seq_type)

    # NTSeqLike will coerce, so this AAs that are all A/T/C/G/U will work here
    # Test that NTSeqLike construction passes when we have guarantees over
    # the type
    if seq_type == "NT":
        s_temp = ntSeqLike(sequence)
        assert str(s_temp) == str(seq)
        assert s_temp._type == seq_type

    # Test that construction of NTSeqLike fails when trying to make an NTSeqLike
    # from a sequence that contains non-NT characters, e.g. AA seqs
    else:
        if not (is_NT(sequence)):
            with pytest.raises(TypeError):
                str(ntSeqLike(sequence)) == str(seq)


@given(string_sequences())
def test_AASeqLike(sequence_type_and_alphabet):
    """Test initialization of SeqLike objects from strings."""
    sequence, seq_type, alphabet = sequence_type_and_alphabet
    seq = SeqLike(sequence, alphabet=alphabet, seq_type=seq_type)

    # AASeqLike will coerce, and every NT sequence is a valid AA sequence
    s_temp = aaSeqLike(sequence)
    assert str(s_temp) == str(seq)
    assert s_temp._type == "AA"
