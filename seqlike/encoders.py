from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder
from .alphabets import NT, AA, STANDARD_AA, STANDARD_NT


def index_encoder_from_alphabet(alphabet):
    """Return a OrdinalEncoder from the tokens in alphabet while preserving order.

    :param alphabet: a iterable of unique tokens
    :returns: OrdinalEncoder
    """
    categories = [[letter for letter in alphabet]]
    fit_list = [[letter] for letter in alphabet]
    return OrdinalEncoder(dtype=float, categories=categories).fit(fit_list)


def onehot_encoder_from_alphabet(alphabet):
    """Return a OneHotEncoder from the tokens in alphabet while preserving order.

    :param alphabet: a iterable of unique tokens
    :returns: OneHotEncoder
    """
    categories = [[letter for letter in alphabet]]
    fit_list = [[letter] for letter in alphabet]
    return OneHotEncoder(dtype=float, sparse=False, categories=categories).fit(fit_list)


ENCODERS = {
    "NT": {
        "index": {
            "standard": index_encoder_from_alphabet(STANDARD_NT),
            "full": index_encoder_from_alphabet(NT),
        },
        "onehot": {
            "standard": onehot_encoder_from_alphabet(STANDARD_NT),
            "full": onehot_encoder_from_alphabet(NT),
        },
    },
    "AA": {
        "index": {
            "standard": index_encoder_from_alphabet(STANDARD_AA),
            "full": index_encoder_from_alphabet(AA),
        },
        "onehot": {
            "standard": onehot_encoder_from_alphabet(STANDARD_AA),
            "full": onehot_encoder_from_alphabet(AA),
        },
    },
}
