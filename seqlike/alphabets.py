"""SeqLike catalog of alphabets."""

import string
from Bio.Data.IUPACData import ambiguous_dna_values, ambiguous_rna_values
from Bio.Data.IUPACData import protein_letters, extended_protein_letters

gap_letter = "-"
stop_letter = "*"
generic_protein_letter = "X"
generic_nt_letter = "N"

# use full alphabet as default for functions that require one
every_letter_alphabet = string.ascii_uppercase

# The rationale for this ordering is that the gap character and standard symbols (4 bases / 20 amino acids) should come first,
# followed by the extra letters.  If we were to use something like alphanumeric ordering, then the standard and full alphabets
# would be mutually incompatible.

STANDARD_NT = gap_letter + "ACGTU" + generic_nt_letter
NT = STANDARD_NT + "BDHKMRSVWY"
STANDARD_AA = stop_letter + gap_letter + protein_letters + generic_protein_letter
AA = STANDARD_AA + "BJOUZ"

STANDARD_NT_SET = set(STANDARD_NT)
NT_SET = set(NT)
STANDARD_AA_SET = set(STANDARD_AA)
AA_SET = set(AA)


# combine ambiguous_dna_values and ambiguous_rna_values into one dict
# :sa: https://stackoverflow.com/questions/1495510/combining-dictionaries-of-lists-in-python
def merge_dicts_of_str(d1, d2, ignore_keys=None):
    if ignore_keys is None:
        ignore_keys = list()
    keys = set(d1).union(d2) - set(ignore_keys)
    return dict((k, "".join(sorted(set(d1.get(k, "") + d2.get(k, ""))))) for k in keys)


ambiguous_nt_values = merge_dicts_of_str(ambiguous_dna_values, ambiguous_rna_values, ignore_keys="X")


# this seq->set->upper->set is necessary to avoid Seq.upper() errors
# (fails for string alphabets by trying to apply alphabet._upper())
# while extracting just the sequence letters (str(SeqRecord) returns
# a string description of the SeqRecord ID, name, etc
def is_NT(sequence):
    # str, Seq, SeqRecord or SeqLike
    return _is_seqtype(sequence, NT_SET)


def is_AA(sequence):
    return _is_seqtype(sequence, AA_SET)


def is_STANDARD_AA(sequence):
    return _is_seqtype(sequence, STANDARD_AA_SET)


def is_STANDARD_NT(sequence):
    return _is_seqtype(sequence, STANDARD_NT_SET)


def _is_seqtype(sequence, seq_letters):

    # seqlike
    if hasattr(sequence, "_seqrecord"):
        sequence = sequence._seqrecord.seq._data

    # seqrecord
    elif hasattr(sequence, "seq"):
        # seqrecord was initialized from a Seq
        try:
            sequence = sequence.seq._data
        # seqrecord was initialized from a string
        except AttributeError:
            sequence = sequence.seq
    # seq
    elif hasattr(sequence, "_data"):
        sequence = sequence._data

    if isinstance(sequence, bytes):
        sequence = sequence.decode()
    sequence = sequence.upper()

    # The meat of the logic lies here.
    return set(sequence).issubset(seq_letters)


def parse_alphabet(alphabet: str) -> str:
    """
    This function parses and validates the 'alphabet' parameter of a SeqLike.

    :param alphabet: str specifying 'NT', 'DNA', 'RNA', or 'AA', case insensitive.
    :returns: either the NT or AA alphabet string.
    """
    # parse string designation to desired alphabet
    if isinstance(alphabet, str):
        alphabet = alphabet.upper()
        assert alphabet in ["NT", "DNA", "RNA", "AA"], "Invalid alphabet!"
    if alphabet in ["DNA", "NT", "RNA"]:
        return NT
    else:
        return AA
