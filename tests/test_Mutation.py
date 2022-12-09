"""Test Mutation class and its features."""
import pytest
from hypothesis import given, strategies as st, assume
from seqlike.alphabets import STANDARD_AA
from seqlike.Mutation import Mutation, Substitution, Deletion, Insertion


@st.composite
def wt_letters(draw):
    """Composite strategy to draw wild-type letters.

    Note here that an empty string is perfectly legitimate,
    as mutations need not necessarily have a wild-type letter specified!
    """
    possibilities = [""] + list(STANDARD_AA)
    return draw(st.sampled_from(possibilities))


@st.composite
def mutant_letters(draw):
    """Composite strategy to draw mutant letters."""
    return draw(st.sampled_from(STANDARD_AA))


@given(wt_letter=wt_letters(), position=st.integers(), mutant_letter=mutant_letters())
def test_Mutation_constructor(wt_letter, position, mutant_letter):
    """Test Mutation (and child class) constructors when passing in wt, pos, mut."""
    assume(wt_letter != "-")
    assume(position > 0)
    position += 1
    mutation = Mutation(f"{wt_letter}{position}{mutant_letter}")
    if wt_letter == "":
        assert mutation.wt_letter is None
    assert mutation.position == position
    assert mutation.mutant_letter == mutant_letter

    string_rep = str(mutation)
    if wt_letter == "":
        assert string_rep == f"{mutation.position}{mutation.mutant_letter}"
    else:
        assert string_rep == f"{mutation.wt_letter}{mutation.position}{mutation.mutant_letter}"

    # Assert equality
    mutation2 = Mutation(f"{wt_letter}{position}{mutant_letter}")
    assert mutation == mutation2


def test_Mutation_magical_constructor():
    """Example-based tests for magical constructors.

    In this test, we are checking that the Mutation constructor
    returns the appropriate child class.
    """

    m1 = Mutation("3R")
    assert isinstance(m1, Substitution)

    m2 = Mutation("^5K")
    assert isinstance(m2, Insertion)

    m3 = Mutation("I4-")
    assert isinstance(m3, Deletion)
