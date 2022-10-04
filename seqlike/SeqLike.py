"""Core SeqLike implementation."""
import itertools

from warnings import warn

import uuid
import warnings
from copy import deepcopy
from functools import reduce
from typing import Callable, Optional, Union

import lazy_loader as lazy

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multipledispatch import dispatch
import pandas as pd

from .alphabets import (
    AA,
    NT,
    STANDARD_AA,
    STANDARD_NT,
    gap_letter,
    is_AA,
    is_NT,
    is_STANDARD_AA,
    is_STANDARD_NT,
)
from .encoders import index_encoder_from_alphabet, onehot_encoder_from_alphabet, array_to_string
from .SequenceLike import SequenceLike
from .utils import (
    add_seqnums_to_letter_annotations,
    ungap,
)
from .Mutation import Mutation, Substitution, Deletion, Insertion
from .MutationSet import MutationSet

np = lazy.load("numpy")


# TODO: Do we want to do some arithmetic on types here?
ArrayType = Union[
    list,
    np.ndarray,
]
StringLikeType = Union[str, Seq, SeqRecord]
SeqLikeType = Union[ArrayType, StringLikeType, "SeqLike"]

seqrecord_attrs = (
    "id",
    "seq",
    "name",
    "description",
    "dbxrefs",
    "features",
    "annotations",
    "letter_annotations",
)

SEQTYPE_TO_TYPE_MAPPING = {"DNA": "NT", "RNA": "NT", "NT": "NT", "AA": "AA"}
NT_TYPE = ["NT", "DNA", "RNA"]
AA_TYPE = ["AA"]


class SeqLike(SequenceLike):
    """
    An omnibus object for various representations of biological sequences.

    This class provides a simple way to interconvert between array,
    string, Seq, and SeqRecord representations of biological sequences.
    In general SeqLike objects try to act something like "dual SeqRecords",
    where we have both nucleotide and amino acid representations.
    All operations tend to act on the current selected form (NT or AA),
    and we strive to keep the representations in sync whenever possible.

    Here's a quick usage example:

    ```python
    from seqlike import SeqLike
    example = 'ATCGATC'

    seq_record = SeqLike(example, seq_type="nt").to_seqrecord()
    seq = SeqLike(example, seq_type="nt").to_seq()
    seq_str = SeqLike(example, seq_type="nt").to_str()
    seq_index = SeqLike(example, seq_type="nt").to_index()
    seq_onehot = SeqLike(example, seq_type="nt").to_onehot()
    ```

    Any of the aformentioned representations can be generated or passed in as an input.
    If using onehot or index encodings, they are of the shape:
    (N x NUM_BASES) and (N), respectively.
    We use symbols from our extended nucleotide and protein alphabets.
    These allow for gaps (`-`), stops (`*`),
    and additional amino acid letters.
    For more info, please see the *_onehot_encoder and *_index_encoder objects,
    and specifically, their categories_ attributes.

    We also allow switching between the NT and AA views of the object.

    ```python
    # conversion to AA
    aa_example = SeqLike(example).aa()

    # conversion to NT
    example_orig = aa_example.nt()
    ```

    The constructor also takes optional keyword arguments
    that are passed on to SeqRecord.
    For the id attribute, if unspecified we generate a random hexadecimal UUID.

    `codon_map` defines back-translation from an AA sequence to a NT sequence.
    This is formatted as a callable, see `codon_tables.py` for more info.

    If the sequence is an NT, then its length must be a multiple of 3 to be translated.
    If it's not, then all AA attributes are None.
    Likewise, if the initializing sequence is an AA,
    but no codon_map was passed in, we can't generate a NT sequence.
    We defer and calculate these on the fly under the hood if possible.
    """

    def __init__(
        self,
        sequence: SeqLikeType,
        seq_type: str,
        alphabet: Optional[str] = None,
        codon_map: Optional[Callable] = None,
        **kwargs,
    ):
        """Initialize a SeqLike object.

        :param sequence: String, Seq, SeqRecord, SeqLike, List, NumPy ndarray, PyTorch Tensor
        :param seq_type: one of 'NT', 'DNA', 'RNA' or 'AA', case insensitive
        :param alphabet: one of 'standard', 'full', or a custom string of alphabet letters
        :param codon_map: a callable, see codon_tables.py for more info.
        :param **kwargs: keyword arguments for internal SeqRecord attribute
        """

        (
            _type,
            _aa_record,
            _nt_record,
            alphabet,
            codon_map,
            _index_encoder,
            _onehot_encoder,
        ) = _construct_seqlike(sequence, seq_type, alphabet, codon_map, **kwargs)

        self._type = _type
        self.alphabet = alphabet
        self._index_encoder = _index_encoder
        self._onehot_encoder = _onehot_encoder
        self._aa_record = _aa_record
        self._nt_record = _nt_record
        self._seqrecord = self._aa_record if _type == "AA" else self._nt_record
        self.codon_map = codon_map

    @property
    def sequence(self):
        return self._seqrecord.seq

    def nt(self, auto_backtranslate=True, **kwargs) -> "SeqLike":
        """
        This method returns the NT view of the SeqLike object.

        The method will automagically back-translate the `._aa_record`
        of its SeqLike object if the `._nt_record` does not exist.

        :param auto_backtranslate: Whether or not
            to automagically back-translate an AA sequence into an NT sequence.
            Defaults to True.
        :param kwargs: Passed through to `.back_translate()`.
            Can include `codon_map` to specify which codon map to use.
        :returns: The NT view of the object.
        """
        # Start with auto-back-translation
        if self._aa_record and self._nt_record is None and auto_backtranslate:
            return self.back_translate(**kwargs)

        if self._type == "NT":
            return deepcopy(self)
        else:
            return swap_representation(self)

    def aa(self, auto_translate=True, **kwargs) -> "SeqLike":
        """
        Return the amino acid view of the SeqLike object.

        The method will automagically translate the `._nt_record`
        of its SeqLike object if the `._aa_record` does not exist.

        We preserve all NT record attributes for convenience in sequence manipulation
        tasks where the starter sequence is an NT and we merely want the AA form.
        The default behaviour, thus, has `id=True`, `name=True`, `description=True`,
        `annotations=True`, and `dbxrefs=True`.

        To change the behaviour,
        you can set any of those to False in `.aa()`, for example:

        ```python
        seq.aa(id=False)
        ```

        :param auto_translate: Whether to automagically translate
            an NT sequence into an AA sequence.
            Defaults to True.
        :param kwargs: These kwargs are passed into BioPython SeqRecord's
            `.translate()` method.  The default behaviour is set to `id=True`,
            `name=True`, `description=True`, `annotations=True`, and `dbxrefs=True`.
        :returns: copy of self with sequence object as amino acid sequence.
        """
        # Start with auto-translation
        if self._nt_record and self._aa_record is None and auto_translate:
            translate_kwargs = dict(id=True, name=True, description=True, annotations=True, dbxrefs=True)
            translate_kwargs.update(kwargs)
            return self.translate(**translate_kwargs)

        # Return based on _type.
        if self._type == "AA":
            return deepcopy(self)
        return swap_representation(self)

    # ------------------------- Methods to inter-convert Sequence to various Biopython formats -------------------------
    # Function to convert the input sequence object to Seq object
    def to_seq(self) -> Seq:
        """
        Convert the SeqLike object (string, Seq, SeqRecord) to a BioPython Seq object.

        :returns: A Seq object.
        """
        return Seq(str(self._seqrecord.seq))

    # Function to convert input sequence to SeqRecord object
    def to_seqrecord(self, **kwargs) -> SeqRecord:
        """
        Convert the SeqLike object (string, Seq, SeqRecord) to a BioPython SeqRecord.

        :param **kwargs: Used to set the attributes of the SeqRecord.
            Overrides the original SeqRecord's attributes.
        :returns: A SeqRecord object.
        """
        rec = self._seqrecord
        newrec = SeqRecord(
            Seq(str(rec.seq)),
            id=rec.id,
            name=rec.name,
            description=rec.description,
            dbxrefs=rec.dbxrefs[:],
            features=rec.features[:],
            annotations=rec.annotations.copy(),
            letter_annotations=rec.letter_annotations.copy(),
        )
        # overwrite newrec values with keyword args
        for key, val in kwargs.items():
            setattr(newrec, key, val)
        return newrec

    # ------------------------- Functions to apply callables -------------------------
    def translate(self, **kwargs) -> "SeqLike":
        """Translate the SeqLike._nt_record and add a ._aa_record.

        Returns the `.aa()` view of the SeqLike object.

        :param kwargs: These kwargs are passed into
            BioPython SeqRecord's `.translate()` method.
            They can be used to set SeqRecord properties.
        :returns: A translated SeqLike object.
        :raises ValueError: When we're trying to translate a SeqLike object
            that doesn't have an NT record.
        :raises TypeError: When we're trying to translate a nucleotide sequence
            that isn't a multiple of 3.
        """
        sc = deepcopy(self)
        if sc._nt_record is None:
            raise ValueError(
                "Oops! It looks like you're trying to translate a SeqLike object "
                "that doesn't have a nucleotide record set. "
                "Unfortunately this would be semantically incorrect. "
                "Please ensure that your SeqLike has a `._nt_record` SeqRecord "
                "before calling on `.translate()`."
            )

        if len(sc) % 3 != 0:
            raise TypeError(
                "Oh no! It looks like you're trying to translate a nucleotide sequence "
                "whose length is not a multiple of 3. "
                "As a safeguard, SeqLike objects do not allow this to happen. "
            )
        sc._nt_record.annotations["molecule_type"] = "DNA"
        sc._aa_record = record_from(sc._nt_record.translate(gap=gap_letter, **kwargs))
        return sc.aa()

    def back_translate(self, codon_map: Callable = None, **kwargs) -> "SeqLike":
        """This method backtranslates the current AA sequence and returns an
        NT sequence.

        We expect that self.codon_map will be defined, and if not, we
        expect a callable and its arguments to be passed into this
        method, as either the first positional argument or by using
        keyword arguments.  Any other arguments for the back
        translator should use keyword arguments.

        If `codon_map` is one of the keyword arguments or the sole
        positional argument, we override self.codon_map, with the
        associated value, and pass all kwargs through to it, along
        with self.  Otherwise we use self.codon_map.

        Finally, we use self.apply to apply the codon_map callable to
        the current sequence.

        :param codon_map: A codon map callable.
        :param **kwargs: Passed through to the codon_map function.
        :returns: A new NT seqlike.
        :raises AttributeError: if no codon map is passed in
            when the SeqLike's `codon_map` is also not set.
        """
        # normally, we wouldn't have to do this sort of thing, we could
        # just use self.__dict__ in the codon_map function.

        # but we might want to use a different codon_map than the original
        # object, so we have to copy and set codon_map in the copy as to not
        # modify in-place.

        codon_map = codon_map or self.codon_map

        if codon_map:
            validate_codon_map(codon_map)
            sc = self.apply(codon_map, **kwargs)  # returns a deepcopied SeqLike

            # TODO: Change this to an if/raise block.
            assert (
                isinstance(sc, SeqLike) and sc._type == "NT"
            ), f"Backtranslating function must return an NT SeqLike, type was {type(sc)}"
        else:
            raise AttributeError(
                "No callable passed and self.codon_map not set!  "
                "Please set the codon_map attribute or pass in a "
                "callable using the `codon_map` keyword."
            )

        return sc.nt()

    def reverse_complement(self, **kwargs) -> "SeqLike":
        """Return reverse complement of NT sequence.

        Record this operation in annotations['reversed'].
        If sequence is currently AA, raise exception like Bio.Seq.Seq.complement()

        :param **kwargs: Not currently used.
        :returns: Reverse-complemented SeqLike object.
        :raises ValueError: when trying to reverse-complement an AA sequence.
            This is a semantically invalid operation.
        """
        if self._type == "AA":
            raise ValueError("Proteins do not have complements!")

        if hasattr(self, "annotations"):
            annotations = self.annotations.copy()
        else:
            annotations = dict()

        if "reversed" in annotations:
            annotations["reversed"] = not annotations["reversed"]
        else:
            annotations["reversed"] = True

        _nt_record = self._nt_record.reverse_complement(
            id=True, name=True, description=True, annotations=annotations, dbxrefs=True
        )

        s = SeqLike(
            _nt_record,
            seq_type=self._type,
            alphabet=self.alphabet,
            codon_map=self.codon_map,
        )
        return s

    # ------------------------- Functions to format sequences -------------------------
    def ungap(self) -> "SeqLike":
        """Remove gap characters and return a new SeqLike.

        :returns: An ungapped SeqLike with no gap characters.
        """
        new_seq = ungap(self.to_seqrecord())

        return SeqLike(
            new_seq,
            seq_type=self._type,
            alphabet=self.alphabet,
            codon_map=self.codon_map,
            index_encoder=self._index_encoder,
            onehot_encoder=self._onehot_encoder,
        )

    def pad_to(self, new_length: int, pad_char: str = gap_letter, mode: str = "right") -> "SeqLike":
        """This method returns a new SeqLike, with gap characters added to the left
        or right to make a sequence of the specified length. Note that this method
        retains only the current view!

        :param new_length: An integer >= than current length
        :param pad_char: The padding character
        :param mode: one of 'left' or 'right'
        :returns: new SeqLike
        """

        assert isinstance(new_length, (int, np.int0, np.int8, np.int16, np.int32, np.int64))
        assert mode in ["left", "right"]

        diff = new_length - len(self)
        assert diff >= 0, (
            f"Current length is {len(self)} and "
            "requested padded length is {new_length}."
            " Can't do negative padding"
        )

        if mode == "left":
            return SeqLike(
                pad_char * diff + self._seqrecord,
                seq_type=self._type,
                alphabet=self.alphabet,
                codon_map=self.codon_map,
                index_encoder=self._index_encoder,
                onehot_encoder=self._onehot_encoder,
            )

        if mode == "right":
            return SeqLike(
                self._seqrecord + pad_char * diff,
                seq_type=self._type,
                alphabet=self.alphabet,
                codon_map=self.codon_map,
                index_encoder=self._index_encoder,
                onehot_encoder=self._onehot_encoder,
            )

    def seq_num_to_idx(self, list_of_seqnums):
        """Convert seqnum to idx.

        TODO: This docstring might need to be made better. tag Andrew Giessel.

        :param list_of_seqnums: Something.
        :returns: Something.
        """

        def ranger(i):
            """TODO: Docstrings need to be added here. tag Andrew Giessel.

            :param i: Something
            :yields: Something
            """
            for _, b in itertools.groupby(enumerate(i), lambda x: x[1] - x[0]):
                b = list(b)
                if len(b) > 1:
                    yield slice(b[0][1], b[-1][1] + 1)
                else:
                    yield b[0][1]

        seqnums_map = dict((seqnum, i) for i, seqnum in enumerate(self.letter_annotations["seqnums"]))
        seqnums = [seqnums_map[seqnum] for seqnum in list_of_seqnums]
        return list(ranger(seqnums))

    # ------------------------- String-like methods -----------------------------------
    def upper(self):
        """Return uppercase sequences and alphabet.

        <!-- #noqa: DAR201 -->
        """
        seq_copy = deepcopy(self)
        seq_copy._seqrecord = seq_copy._seqrecord.upper()
        if seq_copy._nt_record:
            seq_copy._nt_record = seq_copy._nt_record.upper()
        if seq_copy._aa_record:
            seq_copy._aa_record = seq_copy._aa_record.upper()
        seq_copy._type = seq_copy._type
        seq_copy.alphabet = seq_copy.alphabet.upper()
        return seq_copy

    def slice(self, list_of_seqnums):
        """Use seqnums to sub-index.

        TODO: This docstring needs to be much more improved.

        :param list_of_seqnums: Something
        :returns: Something
        """
        return self.__getitem__(self.seq_num_to_idx(list_of_seqnums))

    # ------------------------- Special methods ---------------------------------------
    def __repr__(self) -> str:
        """Return a representation of the sequence.

        If both nucleotide and an associated amino acid sequence are present,
        we show both.

        <!-- #noqa: DAR201 -->
        """
        if self._type == "NT":
            return f"*** NT: {self._nt_record.__repr__()} \n\nAA: {self._aa_record.__repr__()}"
        return f"NT: {self._nt_record.__repr__()} \n\n*** AA: {self._aa_record.__repr__()}"

    def __setattr__(self, name, value):
        """Set attribute value

        We prioritize setting attribute in _seqrecord before the SeqLike object.

        <!-- #noqa: DAR101 -->
        <!-- #noqa: DAR201 -->
        """

        if name in seqrecord_attrs:
            object.__setattr__(self._seqrecord, name, value)
        else:
            object.__setattr__(self, name, value)

    def __getattr__(self, name):
        """Called if the attribute does not already exist in SeqLike

        Please see [the official Python reference][pyref] for more information.

        [pyref]: https://docs.python.org/3/reference/datamodel.html

        :param name: Attribute to return.
        :raises AttributeError: if the attribute is not a SeqLike
            or SeqRecord attribute.

        <!-- #noqa: DAR201 -->
        """
        if name == "__setstate__":  # pertains to pickling
            raise AttributeError

        if name in seqrecord_attrs:
            return getattr(self._seqrecord, name)
        else:
            raise AttributeError("%s not an attribute of SeqLike or SeqLike._seqrecord" % name)

    def __dir__(self):
        """Override for dir() of a SeqLike's attributes.

        :returns: List of SeqLike attributes and SeqRecord attributes.
        """
        return super().__dir__() + list(seqrecord_attrs)

    def __getitem__(self, index) -> "SeqLike":
        """
        __getitem__ implementation.

        Slicing on SeqLike will slice the associated sequence
        and return a sub-sequence SeqLike object.
        We try to cast the SeqLike to NT if possible,
        index,
        and then create a new sequence out of that.
        If the sub-sequence is not a multiple of three in length,
        then there will be no AA sequence associated with the new SeqLike,
        due to ambiguity of coding frame.
        This is simply the default behavior of the SeqLike constructor.
        Note that this *can* result in a new NT with AA that is out of frame,
        the constructor doesn't know or care.

        By convention, the indices here are in the "units" of the
        current type of SeqLike, i.e. in NTs if DNA or RNA, or AAs if AA

        :param index: integer or slice to access parts of the sequence.
        :returns: A new SeqLike of the same type sliced to the index.
        """
        index = [index] if isinstance(index, (int, slice)) else index

        # slicing is a type-preserving
        seqlike_kwargs = dict(
            seq_type=deepcopy(self._type),
            alphabet=deepcopy(self.alphabet),
            codon_map=deepcopy(self.codon_map),
        )

        all_seqlikes = []
        for idx in index:
            if isinstance(idx, int):
                idx = slice(idx, idx + 1, 1)

            # Handle case when both NT and AA are present:
            # - We need to slice BOTH the _nt_record and _aa_records.
            if self._nt_record is not None and self._aa_record is not None:
                # Firstly, determine the primary and secondary records to slice based on _type
                _aa_record = deepcopy(self._aa_record)
                _nt_record = deepcopy(self._nt_record)

                if self._type == "AA":
                    # Slice the AA followed by 3x boundaries.for NT
                    _aa_record = _aa_record[idx]

                    start = idx.start or 0
                    stop = idx.stop or len(self)
                    step = idx.step or 1

                    new_nt_record = None
                    for i in range(start, stop, step):
                        if new_nt_record is None:
                            new_nt_record = _nt_record[i * 3 : (i + 1) * 3]
                        else:
                            new_nt_record = new_nt_record + _nt_record[i * 3 : (i + 1) * 3]

                    sliced = SeqLike(_aa_record, **seqlike_kwargs)
                    sliced._nt_record = new_nt_record

                elif self._type == "NT":
                    _nt_record = deepcopy(self._nt_record)[idx]
                    try:
                        _aa_record = _nt_record.translate()
                    except Exception as e:
                        warnings.warn(e)
                    sliced = SeqLike(_nt_record, **seqlike_kwargs)
                    sliced._aa_record = _aa_record

            # _aa_record is only present
            elif self._type == "AA" and self._nt_record is None:
                sliced = SeqLike(deepcopy(self._aa_record)[idx], **seqlike_kwargs)

            # _nt_record is only present
            elif self._type == "NT" and self._aa_record is None:
                sliced = SeqLike(deepcopy(self._nt_record)[idx], **seqlike_kwargs)

            all_seqlikes.append(sliced)

        # all_seqlikes is now a list of SeqLikes that we need to concatenate
        return reduce((lambda x, y: x + y), all_seqlikes)

    def __add__(self, other):
        """Add sequence to another sequence.

        Mimics behavior of SeqRecord.__add__.

        For magical behaviour, we assume that other is of the same _type as the self._type.
        Doing so allows us to do:

        ```python
        # `a` is SeqLike, `b` is of variable type
        a + b
        ```

        without specifying any further information.

        :param other: A SeqLike type object.
        :returns: The added SeqLike object.
        """

        return _add(self, other)

    def __radd__(self, other: "SeqLike"):
        """Add another sequence or string to this sequence from the left.

        Implementation for:

        ```python
        self = SeqLike(...)
        other = SeqLike(...)
        result = other + self
        ```

        Mimics behavior of SeqRecord.__radd__.

        <!-- #noqa: DAR101 -->
        <!-- #noqa: DAR201 -->
        """
        # for summation of SeqLike using sum()
        if other == 0:
            return self
        if isinstance(other, SeqLike):
            raise RuntimeError("This should have happened via the __add__ of the other SeqLike being added!")

        # Assume it is a string, Seq, or SeqRecord.
        # Note can't transfer any per-letter-annotations
        return SeqLike(
            other + self.to_seqrecord(),
            seq_type=self._type,
            alphabet=self.alphabet,
            codon_map=self.codon_map,
            id=self._seqrecord.id,
            name=self._seqrecord.name,
            description=self._seqrecord.description,
        )

    def __sub__(self, other):
        return _sub(self, other)

    def __deepcopy__(self, memo):
        """Deepcopy implementation.

        <!-- #noqa: DAR101 -->
        <!-- #noqa: DAR201 -->
        """
        seq_copy = SeqLike(
            self.to_seqrecord(),
            seq_type=self._type,
            alphabet=self.alphabet,
            codon_map=deepcopy(self.codon_map),
        )
        seq_copy._nt_record = deepcopy(self._nt_record)
        seq_copy._aa_record = deepcopy(self._aa_record)
        return seq_copy

    def scan(self, mutant_letter: str):
        """Scan a substitution mutation over the sequence."""
        mutants = []
        for i in range(1, len(self) + 1):
            mutants.append(self + Substitution(f"{i}{mutant_letter}"))
        return pd.Series(mutants)


def ntSeqLike(
    sequence: SeqLikeType, alphabet: Optional[str] = None, codon_map: Optional[Callable] = None, **kwargs
) -> SeqLike:
    """
    Entrypoint function for generating a SeqLike of seq_type=='nt', with the same call signature as the SeqLike
    class.  Will coerce the sequence to NT.  Please see SeqLike for more info.
    """
    try:
        if not kwargs["seq_type"].upper() in ["NT", "RNA", "DNA"]:
            warn(
                f"Trying to initialize an NT SeqLike, but seq_type is set to {kwargs['seq_type']}.  Coercing seq_type to NT"
            )
    except KeyError:
        pass
    kwargs["seq_type"] = "NT"
    return SeqLike(sequence, alphabet=alphabet, codon_map=codon_map, **kwargs)


def aaSeqLike(
    sequence: SeqLikeType, alphabet: Optional[str] = None, codon_map: Optional[Callable] = None, **kwargs
) -> SeqLike:
    """
    Entrypoint function for generating a SeqLike of seq_type=='aa', with the same call signature as the SeqLike
    class.  Will coerce the sequence to AA.  Please see SeqLike for more info.
    """
    try:
        if not kwargs["seq_type"].upper() in ["AA"]:
            warn(
                f"Trying to initialize an AA SeqLike, but seq_type is set to {kwargs['seq_type']}.  Coercing seq_type to AA"
            )
    except KeyError:
        pass
    kwargs["seq_type"] = "AA"
    return SeqLike(sequence, alphabet=alphabet, codon_map=codon_map, **kwargs)


@dispatch(SeqLike, (str, type(None)), (str, type(None)), (object, type(None)))
def _construct_seqlike(sequence, seq_type, alphabet, codon_map, **kwargs) -> tuple:
    """Return attributes to set a SeqLike object.

    When we pass in a SeqLike object,
    all relevant attributes are deep-copied over.

    <!-- #noqa: DAR101 -->
    <!-- #noqa: DAR201 -->
    """
    _type = deepcopy(sequence._type)
    _aa_record = record_from(deepcopy(sequence._aa_record), **kwargs)
    _nt_record = record_from(deepcopy(sequence._nt_record), **kwargs)

    if seq_type is None:
        _type = deepcopy(sequence._type)
    if alphabet is None:
        alphabet = deepcopy(sequence.alphabet)
    if codon_map is None:
        codon_map = deepcopy(sequence.codon_map)

    _index_encoder = deepcopy(sequence._index_encoder)
    _onehot_encoder = deepcopy(sequence._onehot_encoder)

    return (
        _type,
        _aa_record,
        _nt_record,
        alphabet,
        codon_map,
        _index_encoder,
        _onehot_encoder,
    )


@dispatch(
    ArrayType.__args__ + StringLikeType.__args__,  # TODO: Figure out a way to include torch tensors w/o requiring torch
    (str),
    (str, type(None)),
    (object, type(None)),
)
def _construct_seqlike(sequence, seq_type, alphabet, codon_map, **kwargs) -> tuple:
    """Return attributes to set a SeqLike object from str, Seq, and SeqRecord objects.

    <!-- #noqa: DAR101 -->
    <!-- #noqa: DAR201 -->
    """
    validate_codon_map(codon_map)

    # Coerce uppercase for `alphabet` and `seq_type`
    alphabet = alphabet.upper() if alphabet is not None else alphabet
    seq_type = seq_type.upper()

    _type, alphabet = determine__type_and_alphabet(seq_type, alphabet, sequence)

    # Get the encoders - both one-hot and index.
    _index_encoder = index_encoder_from_alphabet(alphabet)
    _onehot_encoder = onehot_encoder_from_alphabet(alphabet)

    # Build the _aa_record or _nt_record attribute.
    validate_sequence(sequence, _type)
    seqrecord = record_from(
        sequence,
        _index_encoder=_index_encoder,
        _onehot_encoder=_onehot_encoder,
        **kwargs,
    )

    _aa_record = None if _type == "NT" else seqrecord
    _nt_record = seqrecord if _type == "NT" else None

    return (
        _type,
        _aa_record,
        _nt_record,
        alphabet,
        codon_map,
        _index_encoder,
        _onehot_encoder,
    )


def validate_codon_map(codon_map) -> None:
    """Validate that codon_map is a callable.

    :param codon_map: A codon map callable.
    :raises TypeError: when the codon map is not a Callable.
    """
    err_msg = (
        "An explicitly passed-in codon_map must be a callable, for e.g., "
        "the output of codon_table_to_codon_map(codon_table). "
        "Did you pass in a codon table dictionary instead? "
        "If so, please pass the dictionary through `codon_table_to_codon_map` first "
        "and pass in the resulting function to the `codon_map` argument."
    )
    if codon_map is not None:
        if not callable(codon_map):
            raise TypeError(err_msg)


@dispatch(StringLikeType.__args__, str)
def validate_sequence(sequence, _type) -> None:
    """Validate str, Seq, or SeqRecord sequence objects.

    :raises TypeError: when the sequence is an invalid for its sequence type.

    <!-- #noqa: DAR101 -->
    """
    validation_func = {
        "NT": is_NT,
        "AA": is_AA,
    }

    err_msg = {
        "NT": "Invalid DNA or RNA sequence!",
        "AA": "Invalid protein sequence!",
    }
    if not validation_func[_type](sequence):
        raise TypeError(err_msg[_type])


@dispatch(ArrayType.__args__, str)
def validate_sequence(sequence, _type) -> None:
    """Validate array-like sequence representations.

    :raises ValueError: when the shape of the sequence is not correct.

    <!-- #noqa: DAR101 -->
    """
    sequence = np.asarray(sequence, dtype=float)
    if sequence.ndim not in (1, 2):
        raise ValueError(
            "Numeric representations must be 1d (index) or 2d (onehot). "
            f"However, the shape of the sequence provided is {sequence.shape}."
        )


@dispatch(SeqLike, str)
def validate_sequence(sequence, _type):
    """We do not need to validate SeqLike sequences. (I think.)

    <!-- #noqa: DAR101 -->
    """
    pass


def swap_representation(s: SeqLike) -> SeqLike:
    """Swap representation of a SeqLike object from NT to AA.

    :param s: The SeqLike for which to swap representation.
    :returns: A SeqLike with swapped representation.
    :raises ValueError: When either the _aa_record or _nt_record is missing.
    """
    if s._aa_record is None or s._nt_record is None:
        # Raise an informative error message.
        raise ValueError(
            "Oops! It looks like the SeqLike object is missing "
            "one of `._aa_record` or `._nt_record`."
            "\n\n"
            "Here are the values for you to inspect: \n"
            f"- SeqLike._aa_record: {s._aa_record}\n"
            f"- SeqLike._nt_record: {s._nt_record}"
            "\n\n"
            "Without both representations present, we can't swap views. "
            'If your SeqLike object is of sequence type "NT", '
            "please `.translate()` it first; "
            'alternatively if your SeqLike object is of sequence type "AA", '
            "please `.back_translate()` it first."
        )
    sc = deepcopy(s)  # sc == "seq copy"

    # swap the _type
    _type = "NT" if s._type == "AA" else "AA"

    # copy over the _aa_record and _nt_record objects.
    _aa_record = deepcopy(s._aa_record)
    _nt_record = deepcopy(s._nt_record)

    # swap out the alphabets
    # When swapping representations, standard aa -> standard nt, aa -> nt, etc.
    alphabet_mapping = {
        STANDARD_AA: STANDARD_NT,
        STANDARD_NT: STANDARD_AA,
        NT: AA,
        AA: NT,
    }

    try:
        alphabet = alphabet_mapping.get(s.alphabet)
    except KeyError:
        raise ValueError(
            "Switching between AA and NT views supported only when using "
            "the standard alphabets provided in the SeqLike library. "
            "This ensures that we can swap the alphabets correctly, "
            "and thus also swap the encoders correctly. "
            "If you are using a custom alphabet, "
            "please set the encoders manually."
        )

    # Obtain new encoders
    _index_encoder = index_encoder_from_alphabet(alphabet)
    _onehot_encoder = onehot_encoder_from_alphabet(alphabet)

    # Now set the attributes correctly.
    sc._type = _type
    sc.alphabet = alphabet
    sc._index_encoder = _index_encoder
    sc._onehot_encoder = _onehot_encoder
    sc._aa_record = _aa_record
    sc._nt_record = _nt_record
    sc._seqrecord = sc._aa_record if _type == "AA" else sc._nt_record
    # Codon map doesn't need to be changed so we do not explicitly set it here.
    return sc


@dispatch(str, type(None), StringLikeType.__args__ + ArrayType.__args__)
def determine__type_and_alphabet(seq_type, alphabet, sequence):
    """Determine _type and alphabet when _type is set and alphabet is None.

    <!-- #noqa: DAR101 -->
    <!-- #noqa: DAR201 -->
    """
    _type = SEQTYPE_TO_TYPE_MAPPING.get(seq_type)
    alphabet = determine_alphabet(_type, sequence)
    return _type, alphabet


@dispatch(str, str, StringLikeType.__args__ + ArrayType.__args__)
def determine__type_and_alphabet(seq_type, alphabet, sequence):
    """Determine _type and alphabet when _type is set and alphabet is set.

    <!-- #noqa: DAR101 -->
    <!-- #noqa: DAR201 -->
    """
    _type = SEQTYPE_TO_TYPE_MAPPING.get(seq_type)
    return _type, alphabet


@dispatch(str, StringLikeType.__args__)
def determine__type(alphabet, sequence) -> str:
    """Determine _type when alphabet is set.

    This function uses _only_ the sequence to determine the alphabet.
    We prioritize the amino acid representation where possible.

    <!-- #noqa: DAR101 -->
    <!-- #noqa: DAR201 -->
    """
    if alphabet == AA or alphabet == STANDARD_AA:
        return "AA"
    elif alphabet == NT or alphabet == STANDARD_NT:
        return "NT"


@dispatch(str, StringLikeType.__args__ + (SeqLike,))
def determine_alphabet(_type, sequence) -> str:
    """Determine alphabet from _type and sequence.

    Enables us to use _type _and_ sequence to control the alphabet.

    <!-- #noqa: DAR101 -->
    <!-- #noqa: DAR201 -->
    """
    # _type = determine_seq_type(seq_type, sequence)
    alphabet = {"NT": NT, "AA": AA, "DNA": NT, "RNA": NT}.get(_type)

    # Finally, determine if we need to use the _standard_ versions of the alphabet.
    if _type == "NT" and is_STANDARD_NT(sequence):
        alphabet = STANDARD_NT
    elif _type == "AA" and is_STANDARD_AA(sequence):
        alphabet = STANDARD_AA
    return alphabet


@dispatch(type(None))
def record_from(sequence, **kwargs) -> SeqRecord:
    """Passthrough for None.

    :param sequence: None.
    :param **kwargs: Passed through to SeqRecord constructor.
    :returns: A SeqRecord object.
    """
    return sequence


@dispatch(ArrayType.__args__)
def record_from(sequence, **kwargs) -> SeqRecord:
    """Construct SeqRecord from array-like sequences.

    :param sequence: An array-like object.
    :param **kwargs: Passed through to SeqRecord constructor.
    :returns: A SeqRecord object.
    """
    _index_encoder = kwargs.pop("_index_encoder")
    _onehot_encoder = kwargs.pop("_onehot_encoder")
    s: str = array_to_string(sequence, _index_encoder, _onehot_encoder)
    return record_from(s, **kwargs)


@dispatch(str)
def record_from(sequence, **kwargs) -> SeqRecord:
    """Construct SeqRecord from string-like sequences.

    :param sequence: A string object.
    :param **kwargs: Passed through to SeqRecord constructor.
    :returns: A SeqRecord object.
    """
    s: Seq = Seq(sequence)
    return record_from(s, **kwargs)


@dispatch(Seq)
def record_from(sequence, **kwargs) -> SeqRecord:
    """Construct SeqRecord from Seq sequences.

    :param sequence: A Seq object.
    :param **kwargs: Passed through to SeqRecord constructor.
    :returns: A SeqRecord object.
    """
    s: SeqRecord = SeqRecord(sequence, id=str(uuid.uuid4()))
    return record_from(s, **kwargs)


@dispatch(SeqRecord)
def record_from(sequence, **kwargs) -> SeqRecord:
    """Construct SeqRecord from SeqRecord sequences.

    :param sequence: A SeqRecord object.
    :param **kwargs: Passed through to SeqRecord constructor.
    :returns: A SeqRecord object.

    """
    s: SeqRecord = deepcopy(sequence)
    for k, v in kwargs.items():
        setattr(s, k, v)
    s = add_seqnums_to_letter_annotations(s)
    return s


@dispatch(SeqLike)
def record_from(sequence, **kwargs) -> SeqRecord:
    """Construct SeqRecord from SeqLike objects.

    NOTE: This is one that is a bit weird but we might find it still useful.

    :param sequence: A SeqLike object.
    :param **kwargs: Passed through to SeqRecord constructor.
    :returns: A SeqRecord object.
    """
    return record_from(deepcopy(sequence._seqrecord), **kwargs)


@dispatch(SeqLike, (str, Seq, SeqRecord, SeqLike))
def _add(obj: SeqLike, other) -> SeqLike:
    """Add a SeqLike-ish object to a SeqLike object.

    :param obj: The SeqLike object.
    :param other: The SeqLike-ish object.
    :returns: A SeqLike object.
    """
    if not isinstance(other, SeqLike):
        other = SeqLike(
            other,
            seq_type=deepcopy(obj._type),
            alphabet=deepcopy(obj.alphabet),
            codon_map=deepcopy(obj.codon_map),
        )
    if not obj._type == other._type:
        raise ValueError(f"Incompatible sequence types! obj is {obj._type}, other is {other._type}")

    return SeqLike(
        obj.to_seqrecord() + other.to_seqrecord(),
        seq_type=obj._type,
        alphabet=obj.alphabet,
        codon_map=obj.codon_map,
        id=obj._seqrecord.id,
        name=obj._seqrecord.name,
        description=obj._seqrecord.description,
    )


@dispatch(SeqLike, (Substitution, Deletion))
def _add(obj: SeqLike, other: Union[Substitution, Deletion]):
    """Add a Substitution or Deletion object to a SeqLike."""
    other = deepcopy(other)
    other._one_indexed = True
    validate_mutation_position(obj, other)
    new = ""
    for i, l in enumerate(obj.to_str()):
        if i != other.position - other._one_indexed:
            new += l
        elif i == other.position - other._one_indexed:
            new += other.mutant_letter
    return SeqLike(
        new,
        seq_type=deepcopy(obj._type),
        alphabet=deepcopy(obj.alphabet),
        codon_map=deepcopy(obj.codon_map),
        id=deepcopy(obj._seqrecord.id),
        name=deepcopy(obj._seqrecord.name),
        description=deepcopy(obj._seqrecord.description),
    )


@dispatch(SeqLike, Insertion)
def _add(obj: SeqLike, other: Insertion):
    """Add an Insertion object to a SeqLike."""
    other = deepcopy(other)
    other._one_indexed = True
    validate_mutation_position(obj, other)
    new = ""
    for i, l in enumerate(obj.to_str()):
        if i == other.position - other._one_indexed:
            new += other.mutant_letter
        new += l
    # Final check for mutations after the final letter.
    if i + 1 == other.position - other._one_indexed:
        new += other.mutant_letter

    return SeqLike(
        new,
        seq_type=deepcopy(obj._type),
        alphabet=deepcopy(obj.alphabet),
        codon_map=deepcopy(obj.codon_map),
        id=deepcopy(obj._seqrecord.id),
        name=deepcopy(obj._seqrecord.name),
        description=deepcopy(obj._seqrecord.description),
    )


@dispatch(SeqLike, MutationSet)
def _add(obj: SeqLike, other: MutationSet):
    """Add a MutationSet object to a SeqLike."""
    mutations = deepcopy(other.mutations)
    while mutations:
        mutation = mutations.pop(0)
        obj = obj + mutation
        if isinstance(mutation, Insertion):
            mutations = [m + 1 for m in mutations]
    return obj


@dispatch(SeqLike, SeqLike)
def _sub(reference: SeqLike, other: SeqLike) -> MutationSet:
    """Subtract a SeqLike from another SeqLike to obtain a MutationSet diff."""
    reference = deepcopy(reference)
    other = deepcopy(other)
    mutations = []
    seqrecs = pd.Series([reference, other]).seq.align()
    reference = seqrecs.iloc[0]
    other = seqrecs.iloc[1]
    i = 1
    for wt, mut in zip(str(reference), str(other)):
        if wt == mut:
            i += 1
            continue
        if wt == "-":
            mutation = Mutation(f"^{i}{mut}")
            mutations.append(mutation)
            continue
        else:
            mutation = Mutation(f"{i}{mut}")
            mutations.append(mutation)
            i += 1
            continue
    return MutationSet(mutations)


def validate_mutation_position(obj: SeqLike, other: Mutation):
    if other.position > len(obj):
        raise ValueError(
            f"Mutation {other} cannot be applied to the SeqLike object because it is at a position greater than the length of the SeqLike object (length = {len(obj)})!"
        )
    actual_wt = str(obj[other.position - other._one_indexed])
    if other.wt_letter and other.wt_letter != "^" and actual_wt != other.wt_letter:
        raise ValueError(
            f"Mutation {other} is expected to be applied to a SeqLike with WT letter {other.wt_letter} at position {other.position} but found actual letter {actual_wt}!"
        )
