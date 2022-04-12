"""
SequenceLike: a more general class that relaxes the assumption
that the symbols in the sequence are amino or nucleic acids,
but preserves the ability to convert to/from index and onehot encodings,
given a sequence and an alphabet.
"""
from collections.abc import Sequence
import numpy as np
from copy import deepcopy

from .encoders import index_encoder_from_alphabet, onehot_encoder_from_alphabet, array_to_symbols


class SequenceLike(Sequence):
    def __init__(self, sequence, alphabet=None, encoding=None):
        if alphabet is None:
            alphabet = set([x for x in sequence])
        alphabet = sorted(alphabet)

        # Get the encoders - both one-hot and index.
        _index_encoder = index_encoder_from_alphabet(alphabet)
        _onehot_encoder = onehot_encoder_from_alphabet(alphabet)

        if encoding in ["onehot", "index"]:
            sequence = array_to_symbols(sequence, _index_encoder, _onehot_encoder)

        # Set properties all in one block for readability.
        self.alphabet = alphabet
        self._index_encoder = _index_encoder
        self._onehot_encoder = _onehot_encoder
        self.sequence = sequence

    def to_str(self) -> str:
        """
        Convert the SequenceLike object to a string.

        :returns: A string.
        """
        return "".join(self.sequence)

    # ------------------------- Functions to convert a sequence to numerical formats -------------------------
    def to_index(self, dtype=int, encoder=None) -> np.ndarray:
        """
        Convert the SequenceLike object into a index-encoded array.

        Underneath the hood we use the `._index_encoder` to handle transformations,
        unless the `encoder` argument is specified.

        :param dtype: Array data type. Defaults to float.
            Can be changed to any dtype that NumPy arrays allow.
        :param encoder: A sklearn encoder
            that transforms the sequence into an integer-encoded array.
        :returns: index encoded numpy array for the given sequence
        """
        seq_as_array = [[x] for x in self]

        if encoder is not None:
            return encoder.transform(seq_as_array).squeeze().astype(dtype)
        return self._index_encoder.transform(seq_as_array).squeeze().astype(dtype)

    def to_onehot(self, dtype=int, encoder=None) -> np.ndarray:
        """
        Convert the SequenceLike object into a one-hot encoded array.

        Underneath the hood we use the `._onehot_encoder` to handle transformations.

        :param dtype: Array data type. Defaults to float.
            Can be changed to any dtype that NumPy arrays allow.
        :param encoder: A sklearn encoder
            that transforms the sequence into an onehot-encoded array.
        :returns: onehot encoded numpy array for the given sequence.
        """
        seq_as_array = [[x] for x in self]

        if encoder is not None:
            return encoder.transform(seq_as_array).squeeze().astype(dtype)
        return self._onehot_encoder.transform(seq_as_array).squeeze().astype(dtype)

    def apply(self, func, **kwargs):
        """This method applies func to the object and returns a copy, using
        keyword arguments.  Note that we *only* support kwargs, not
        positional arguments.  The callable function must take a
        SeqLike as the parameter `seq`, and will often (but not
        always) return a new SeqLike object.

        While one could simply call `func(seq)`, This enables a fluent
        interface pattern, which lets you chain methods together.
        Subclasses of SequenceLike can/should use this method and other
        pre-built functions to rapidly add functionality by defining a
        method that simply calls `self.apply(self, func)`.

        :param func: The function to apply to the SequenceLike object.
        :param **kwargs: Passed through to `func()`.
        :returns: The result of calling `func()`.
        """
        return func(seq=deepcopy(self), **kwargs)

    def find(self, sub, start=None, end=None) -> int:
        """Return index of first occurrence of sub-sequence.  Must be in this
        SequenceLike's alphabet.

        :param sub: A sequence of valid symbols (i.e. must be found in self.alphabet)
        :param start: optional start argument, follows slice notation
        :param end: optional end argument, follows slice notation
        :returns: Integer location of the substring.
        """
        assert set(sub).issubset(self.alphabet)

        if start is None:
            start = 0
        if end is None:
            end = len(self) - len(sub)

        for i in range(start, end + 1):
            if "".join(self[i : i + len(sub)]) == "".join(sub):
                return i
        return -1

    def count(self, sub, start=None, end=None) -> int:
        """Return count of sub-sequence.  Must be in this
        SequenceLike's alphabet.

        :param sub: A sequence of valid symbols (i.e. must be found in self.alphabet)
        :param start: Argument passed on to the underlying `_seqrecord.Seq` object.
        :param end: Argument passed on to the underlying `_seqrecord.Seq` object.
        :returns: Integer location of the substring.
        """
        if not set(sub).issubset(self.alphabet):
            raise ValueError(f"Subsequence set {set(sub)} is not part of the alphabet {self.alphabet}")

        if start is None:
            start = 0
        if end is None:
            end = len(self) - len(sub)

        count = 0
        for i in range(start, end + 1):
            if "".join(self[i : i + len(sub)]) == "".join(sub):
                count += 1
        return count

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __contains__(self, x: object) -> bool:
        return x in self.sequence

    def __iter__(self):
        return iter(self.sequence)

    def __str__(self) -> str:
        """Cast to a string

        <!-- #noqa: DAR201 -->
        """
        return self.to_str()

    def __deepcopy__(self, memo):
        """Deepcopy implementation.

        <!-- #noqa: DAR101 -->
        <!-- #noqa: DAR201 -->
        """
        return SequenceLike(self.sequence, self.alphabet)

    def __repr__(self) -> str:
        return self.sequence.__repr__()
