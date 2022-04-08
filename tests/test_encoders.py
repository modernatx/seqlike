import sklearn
import numpy as np


from seqlike.encoders import index_encoder_from_alphabet, onehot_encoder_from_alphabet


# def test_ENCODERS():
#     for k, v in ENCODERS.items():
#         assert k in ["AA", "NT"]
#         assert isinstance(v, dict)

#     for seq_type in ["AA", "NT"]:
#         for k, v in ENCODERS[seq_type].items():
#             assert k in ["onehot", "index"]
#             assert isinstance(v, dict)

#             for kk, vv in v.items():
#                 assert kk in ["full", "standard"]
#                 assert isinstance(vv, sklearn.preprocessing._encoders._BaseEncoder)


def test_index_encoder():
    enc = index_encoder_from_alphabet("-ACGTUN")
    assert isinstance(enc, sklearn.preprocessing._encoders.OrdinalEncoder)
    assert np.all(enc.categories_[0] == np.array(["-", "A", "C", "G", "T", "U", "N"]))


def test_onehot_encoder():
    enc = onehot_encoder_from_alphabet("-ACGTUN")
    assert isinstance(enc, sklearn.preprocessing._encoders.OneHotEncoder)
    assert np.all(enc.categories_[0] == np.array(["-", "A", "C", "G", "T", "U", "N"]))
