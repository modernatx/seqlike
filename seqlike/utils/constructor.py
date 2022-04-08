# """Functions that are used in the SeqLike constructor.

# They live outside of the SeqLike.py module to keep it unpolluted.

# All of the functions here ought to be imported into its sibling __init__.py file
# so that it's available at the top of the `seqlike.utils` namespace.
# """
# from copy import deepcopy
# from seqlike.utils.validation import validate_seq_type
# from seqlike.encoders import index_encoder_from_alphabet, onehot_encoder_from_alphabet
# from typing import Any, Union

# # import numpy as np
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from multipledispatch import dispatch
# from seqlike.alphabets import (
#     AA,
#     NT,
#     STANDARD_AA,
#     STANDARD_NT,
#     is_AA,
#     is_NT,
#     is_STANDARD_AA,
#     is_STANDARD_NT,
# )
# from seqlike.utils.sequences import add_seqnums_to_letter_annotations
