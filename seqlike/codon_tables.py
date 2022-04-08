"""
A suite of codon tables to be used.  We've included tables from
python_codon_tables below for reference.
"""

from copy import deepcopy
import numpy as np
from seqlike.SeqLike import SeqLike
from typing import Callable
from python_codon_tables import get_codons_table

CODON_TABLE = {
    "AAA": "K",
    "AAC": "N",
    "AAG": "K",
    "AAT": "N",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AGA": "R",
    "AGC": "S",
    "AGG": "R",
    "AGT": "S",
    "ATA": "I",
    "ATC": "I",
    "ATG": "M",
    "ATT": "I",
    "CAA": "Q",
    "CAC": "H",
    "CAG": "Q",
    "CAT": "H",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "GAA": "E",
    "GAC": "D",
    "GAG": "E",
    "GAT": "D",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "TAA": "*",
    "TAC": "Y",
    "TAG": "*",
    "TAT": "Y",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TGA": "*",
    "TGC": "C",
    "TGG": "W",
    "TGT": "C",
    "TTA": "L",
    "TTC": "F",
    "TTG": "L",
    "TTT": "F",
    "NNN": "X",
    "---": "-",
}

# https://github.com/Edinburgh-Genome-Foundry/codon-usage-tables/blob/master/codon_usage_data/tables/h_sapiens_9606.csv
human_codon_table = get_codons_table("h_sapiens_9606")

# https://github.com/Edinburgh-Genome-Foundry/codon-usage-tables/blob/master/codon_usage_data/tables/s_cerevisiae_4932.csv
yeast_codon_table = get_codons_table("s_cerevisiae_4932")

# https://github.com/Edinburgh-Genome-Foundry/codon-usage-tables/blob/master/codon_usage_data/tables/e_coli_316407.csv
ecoli_codon_table = get_codons_table("e_coli_316407")

random_codon_table = {
    "*": {"TAA": 0.33, "TAG": 0.33, "TGA": 0.33},
    "A": {"GCA": 0.25, "GCC": 0.25, "GCG": 0.25, "GCT": 0.25},
    "C": {"TGC": 0.50, "TGT": 0.50},
    "D": {"GAC": 0.50, "GAT": 0.50},
    "E": {"GAA": 0.50, "GAG": 0.50},
    "F": {"TTC": 0.50, "TTT": 0.50},
    "G": {"GGA": 0.25, "GGC": 0.25, "GGG": 0.25, "GGT": 0.25},
    "H": {"CAC": 0.50, "CAT": 0.50},
    "I": {"ATA": 0.07, "ATC": 0.42, "ATT": 0.51},
    "K": {"AAA": 0.50, "AAG": 0.50},
    "L": {"CTA": 0.166, "CTC": 0.166, "CTG": 0.166, "CTT": 0.166, "TTA": 0.166, "TTG": 0.166},
    "M": {"ATG": 1.0},
    "N": {"AAC": 0.50, "AAT": 0.50},
    "P": {"CCA": 0.25, "CCC": 0.25, "CCG": 0.25, "CCT": 0.25},
    "Q": {"CAA": 0.50, "CAG": 0.50},
    "R": {"AGA": 0.166, "AGG": 0.166, "CGA": 0.166, "CGC": 0.166, "CGG": 0.166, "CGT": 0.166},
    "S": {"AGC": 0.166, "AGT": 0.166, "TCA": 0.166, "TCC": 0.166, "TCG": 0.166, "TCT": 0.166},
    "T": {"ACA": 0.25, "ACC": 0.25, "ACG": 0.25, "ACT": 0.25},
    "V": {"GTA": 0.25, "GTC": 0.25, "GTG": 0.25, "GTT": 0.25},
    "W": {"TGG": 1.0},
    "Y": {"TAC": 0.50, "TAT": 0.50},
}


def codon_table_to_codon_map(codon_table: dict, deterministic: bool = True) -> Callable[[SeqLike], SeqLike]:
    """This is a convenience function takes a codon map as defined by the
    dictionaries above and returns a backtranslation callable that
    takes an AA SeqLike and returns a NT SeqLike.

    By default, the SeqLike back_translate() method takes kwargs.  In
    this simple case, we do not.  The codon_table and deterministic
    flag are baked into the returned callable and can't be changed
    after definition.

    :param codon_table: A nested dictionary that has keys being AAs
        and values codon-probabilty pairs
    :param deterministic: if True we return the first codon and disregard the rest.
    :returns: A callable that takes and returns a SeqLike
    """

    def backtranslator(seq):
        if seq._type != "AA":
            raise TypeError("Sequence must be an AA SeqLike!")
        seq_str = seq.to_str()

        nt = ""
        for aa in seq_str:
            codons, probs = zip(*codon_table[aa].items())

            # we normalize the probabilities
            # most tables are near 1.0, but issues with precision exist
            sum_prob = sum(probs)
            probs = [p / sum_prob for p in probs]

            if deterministic:
                nt += codons[0]
            else:
                nt += np.random.choice(codons, p=probs)

        new_seqlike = SeqLike(
            nt,
            id=seq.id,
            name=seq.name,
            description=seq.description,
            annotations=seq.annotations,
            dbxrefs=seq.dbxrefs,
            seq_type="dna",
            codon_map=seq.codon_map,
        )
        new_seqlike._aa_record = deepcopy(seq._aa_record)
        return new_seqlike

    return backtranslator


# add in support for pseudo-translation
for table in [human_codon_table, yeast_codon_table, ecoli_codon_table, random_codon_table]:
    table["-"] = {"---": 1}
    table["X"] = {"NNN": 1}

human_codon_map = codon_table_to_codon_map(human_codon_table)
yeast_codon_map = codon_table_to_codon_map(yeast_codon_table)
ecoli_codon_map = codon_table_to_codon_map(ecoli_codon_table)
