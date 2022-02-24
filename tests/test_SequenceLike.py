"""Tests for SequenceLike."""
from seqlike.SequenceLike import SequenceLike
from collections import Counter
import pytest


def test_inferred_alphabet():
    """Test that the inferred alphabet is correct when constructing a sequence without an alphabet."""
    sequence = ["AD", "AD", "AD", "ASD", "BC", "BCD"]
    s = SequenceLike(sequence)
    assert s.alphabet == sorted(set(sequence))


def test_SequenceLike():
    """Example-based test for SequenceLike object."""
    sequence = ["AD", "AD", "AD", "BC", "ASD", "BCD"]
    s = SequenceLike(sequence)
    assert s.to_str() == "".join(s for s in sequence)
    assert len(s.to_index()) == len(sequence)
    assert len(s.to_onehot()) == len(sequence)
    assert str(s) == s.to_str()

    # Test that `.count()` works.
    counts = Counter(sequence)
    for k, v in counts.items():
        assert s.count([k]) == v

    with pytest.raises(ValueError):
        s.count("AD")
