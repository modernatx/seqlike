from .alphabets import NT, AA, STANDARD_AA, STANDARD_NT
import numpy as np
from typing import Union


def index_encoder_from_alphabet(alphabet):
    """Return a OrdinalEncoder from the tokens in alphabet while preserving order.

    :param alphabet: a iterable of unique tokens
    :returns: OrdinalEncoder
    """
    from sklearn.preprocessing import OrdinalEncoder

    categories = [[letter for letter in alphabet]]
    fit_list = [[letter] for letter in alphabet]
    return OrdinalEncoder(dtype=float, categories=categories).fit(fit_list)


def onehot_encoder_from_alphabet(alphabet):
    """Return a OneHotEncoder from the tokens in alphabet while preserving order.

    :param alphabet: a iterable of unique tokens
    :returns: OneHotEncoder
    """
    from sklearn.preprocessing import OneHotEncoder

    categories = [[letter for letter in alphabet]]
    fit_list = [[letter] for letter in alphabet]
    return OneHotEncoder(dtype=float, sparse_output=False, categories=categories).fit(fit_list)


def array_to_symbols(sequence: Union[list, np.ndarray], _index_encoder, _onehot_encoder) -> str:
    """Convert array-like sequence representations to a string.

    :raises IndexError: if the alphabet of the sequence is a superset
        of the index encoder and one-hot encoder object.

    <!-- #noqa: DAR101 -->
    <!-- #noqa: DAR201 -->
    """
    sequence = np.asarray(sequence, dtype=float)
    if sequence.ndim == 1:
        try:
            sequence = _index_encoder.inverse_transform(sequence.reshape(-1, 1)).flatten()
        except IndexError:
            raise IndexError(
                "The encoder encountered a bad encoding value. "
                "Ensure that you're using an alphabet "
                "which contains all needed symbols."
            )
    elif sequence.ndim == 2:
        sequence = _onehot_encoder.inverse_transform(sequence).flatten()

    # NOTE: We do not need to check for other dim sizes
    # because we assume that validate_sequence will take care of it.
    return sequence


def array_to_string(sequence: Union[list, np.ndarray], _index_encoder, _onehot_encoder) -> str:
    """Convert array-like sequence representations to a string.

    :raises IndexError: if the alphabet of the sequence is a superset
        of the index encoder and one-hot encoder object.

    <!-- #noqa: DAR101 -->
    <!-- #noqa: DAR201 -->
    """
    return "".join(array_to_symbols(sequence, _index_encoder, _onehot_encoder))


# Defer execution till when it's needed.
# ENCODERS = {
#     "NT": {
#         "index": {
#             "standard": index_encoder_from_alphabet(STANDARD_NT),
#             "full": index_encoder_from_alphabet(NT),
#         },
#         "onehot": {
#             "standard": onehot_encoder_from_alphabet(STANDARD_NT),
#             "full": onehot_encoder_from_alphabet(NT),
#         },
#     },
#     "AA": {
#         "index": {
#             "standard": index_encoder_from_alphabet(STANDARD_AA),
#             "full": index_encoder_from_alphabet(AA),
#         },
#         "onehot": {
#             "standard": onehot_encoder_from_alphabet(STANDARD_AA),
#             "full": onehot_encoder_from_alphabet(AA),
#         },
#     },
# }
