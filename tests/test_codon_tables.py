import pytest

from seqlike import aaSeqLike
from seqlike.codon_tables import codon_table_to_codon_map, human_codon_table


@pytest.mark.parametrize("letter", list("ACDEFGHIKLMNPQRSTVWY"))
def test_codon_table_to_codon_map(letter):
    human_codon_map = codon_table_to_codon_map(human_codon_table)
    codon = aaSeqLike(letter).back_translate(codon_map=human_codon_map)
    assert sorted(human_codon_table[letter].items(), key=lambda x: x[1], reverse=True)[0][0] == str(codon)
