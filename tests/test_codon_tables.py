import pytest

from python_codon_tables import get_codons_table
from seqlike import aaSeqLike
from seqlike.codon_tables import (
    codon_table_to_codon_map,
    human_codon_table,
    ecoli_codon_table,
    sort_codon_table_by_frequency,
)


@pytest.mark.parametrize("codon_table", [human_codon_table, ecoli_codon_table])
def test_sort_codon_table_by_frequency(codon_table):
    sorted_codon_table = sort_codon_table_by_frequency(codon_table)
    for letter, codon_frequencies in codon_table.items():
        assert sorted_codon_table[letter] == dict(sorted(codon_frequencies.items(), key=lambda x: x[1], reverse=True))


@pytest.mark.parametrize("letter", list("ACDEFGHIKLMNPQRSTVWY"))
def test_codon_table_to_codon_map(letter):
    # sampled
    human_codon_map = codon_table_to_codon_map(human_codon_table, deterministic=False)
    codon = aaSeqLike(letter).back_translate(codon_map=human_codon_map)
    assert str(codon) in human_codon_table[letter]

    # deterministic
    human_codon_map = codon_table_to_codon_map(human_codon_table)
    codon = aaSeqLike(letter).back_translate(codon_map=human_codon_map)
    assert list(human_codon_table[letter])[0] == str(codon)

    # deterministic sorted
    human_codon_map = codon_table_to_codon_map(sort_codon_table_by_frequency(human_codon_table))
    codon = aaSeqLike(letter).back_translate(codon_map=human_codon_map)
    assert sorted(human_codon_table[letter].items(), key=lambda x: x[1], reverse=True)[0][0] == str(codon)


@pytest.mark.parametrize("codon_table_name", ["h_sapiens_9606", "s_cerevisiae_4932"])
def test_codon_table(codon_table_name):
    assert ''.join(sorted(human_codon_table.keys())) == '*-ACDEFGHIKLMNPQRSTVWXY'

    codons_table = get_codons_table(codon_table_name)
    assert ''.join(sorted(codons_table.keys())) == '*ACDEFGHIKLMNPQRSTVWY'