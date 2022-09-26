"""Test Mutation class and its features."""
import pytest
from hypothesis import given, strategies as st
from seqlike.alphabets import STANDARD_AA
from seqlike.Mutation import Mutation


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
    mutation = Mutation(wt_letter=wt_letter, position=position, mutant_letter=mutant_letter)
    assert mutation.wt_letter == wt_letter
    assert mutation.position == position
    assert mutation.mutant_letter == mutant_letter

    string_rep = str(mutation)
    assert string_rep == f"{mutation.wt_letter}{mutation.position}{mutation.mutant_letter}"

    # Assert equality
    mutation2 = Mutation(wt_letter=wt_letter, position=position, mutant_letter=mutant_letter)
    assert mutation == mutation2
