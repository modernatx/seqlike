import io
from typing import Optional
import warnings

import numpy as np

# import pandas as pd
import lazy_loader as lazy
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data.IUPACData import extended_protein_values
from PIL import Image


from .alignment_utils import align
from .alphabets import (
    AA,
    NT,
    STANDARD_AA,
    STANDARD_NT,
    ambiguous_nt_values,
    gap_letter,
    generic_nt_letter,
    generic_protein_letter,
    stop_letter,
)
from .draw_utils import draw_alignment, view_alignment, aa_chemistry_simple, nt_simple
from .encoders import onehot_encoder_from_alphabet
from .SeqLike import SeqLike

pd = lazy.load("pandas")


@pd.api.extensions.register_series_accessor("seq")
class SeqLikeAccessor:
    """This class extends the Pandas Series class to include a 'seq'
    namespace, which exposes a number of methods that work on a series
    of SeqLikes.  Because of the decorator, when SeqLike is imported,
    this extension is automatically is in effect. Methods are accessed
    by calling `df["column_name"].seq.method()`

    By convention, we expects that the series is comprised of SeqLike
    objects all of the same type.  They do not need to be the same
    length.

    Most methods will return new Pandas Series with new SeqLikes (as
    copies).
    """

    def __init__(self, pandas_obj):
        self._validate(pandas_obj)
        self._obj = self._match_alphabets(pandas_obj)

    @staticmethod
    def _validate(obj):
        """Static method to make sure we have a 'seqs' column and they're all
        SeqLikes of the same type

        :param obj: A SeqLike object.
        :raises ValueError: When seqlikes do not contain the same `_type` attribute.
        """

        # We assume that the 'seqs' column is all SeqLikes
        is_seqlikes = obj.apply(lambda x: isinstance(x, SeqLike))
        types = obj.drop_duplicates().apply(lambda x: type(x))
        if not all(is_seqlikes):
            raise ValueError(f"Series must contain all SeqLike objects, instead found {types}")

        # All SeqLikes must be of the same type
        if len(set(obj.apply(lambda s: s._type.upper()))) != 1:
            raise ValueError("All SeqLikes must be of the same _type (i.e. 'aa', or 'nt').")

    @staticmethod
    def _match_alphabets(obj):
        """If the alphabets are not identical, set them all to the 'full' alphabet.

        :note: This will clobber any custom alphabets.
        :param obj: The pandas Series to apply this function to.
        :returns: A modified version of the pandas Series
            where each SeqLike object's alphabet has been made "full".
        """

        alphabets = obj.apply(lambda s: s.alphabet)
        if len(alphabets.unique()) != 1:
            # TODO: we might want to convert only those SeqLikes
            # that don't have the same alphabet as the original.
            # Doing so could avoid a bunch of deepcopying.
            warnings.warn(
                "It appears that the sequences here have multiple alphabets. "
                "We are replacing alphabets with the full version (AA/NT) "
                "for the full collection. "
            )
            if obj.iloc[0]._type == "AA":
                return obj.apply(lambda s: SeqLike(s, seq_type=s._type, alphabet=AA))
            else:
                return obj.apply(lambda s: SeqLike(s, seq_type=s._type, alphabet=NT))
        return obj

    def write(self, *args, **kwargs):
        """Simple wrapper on `SeqIO.write`.

        :param *args: Passed into `SeqIO.write`.
        :param **kwargs: Passed on to `SeqIO.write`.
        """
        SeqIO.write([seq.to_seqrecord() for seq in self._obj], *args, **kwargs)

    def plot(self, use_bokeh=True, colorscheme=None, x_scale=1, y_scale=1, *args, **kwargs):
        """Plot the SeqLikes as a multiple sequence alignment.

        All *args and **kwargs parameters mirror
        .draw_utils.draw_alignment and .draw_utils.view_alignment.
        We use .as_alignment() for convenience.

        :param use_bokeh: bool; if True (default), use Bokeh backend if available, otherwise use draw_alignment
        :param colorscheme: ColorScheme, WebLogo based mapping of symbol to color
        :param x_scale: float, x scaling factor used with draw_alignment backend
        :param y_scale: float, y scaling factor used with draw_alignment backend
        :param *args: Passed into `seqlike.alignment_utils.view_alignment`.
        :param **kwargs: Passed into `seqlike.alignment_utils.view_alignment`.
        :returns: a PIL Image object or Bokeh object.
        """
        if colorscheme is None and self._type == "NT":
            colorscheme = nt_simple
        elif colorscheme is None and self._type == "AA":
            colorscheme = aa_chemistry_simple

        if use_bokeh:
            try:
                return view_alignment(self.as_alignment(), colorscheme=colorscheme, *args, **kwargs)
            except:
                pass

        x = draw_alignment(self.as_alignment(), colorscheme=colorscheme, *args, **kwargs)
        return x.resize(size=(int(x.size[0] * x_scale), int(x.size[1] * y_scale)))

    def weblogo(
        self,
        seqnum_labels=None,
        ref_id=None,
        cols=50,
        color_scheme=aa_chemistry_simple,
        logo_font="ArialMT",
        logo_format="png",
        resolution=200,
        **kwargs,
    ):
        """Draw weblogo from sequence alignment with optional labeling for consensus mutations
        :param seqnum_labels: label the weblogo letters with these seqnums
        :param ref_id: derive the weblogo letter labels using seqnums from this reference record
        :param cols: weblogo column width (number of letters)
        :param color_scheme: weblogo color scheme
        :param logo_font: weblogo font
        :param logo_format: weblogo image format ('png', 'eps', 'jpg', 'eps', etc)
        :param resolution: weblogo resolution in DPI
        :param **kwargs: additional weblogo arguments
        :returns: PIL Image object
        """
        import weblogo as wl

        def highlight_consensus(labels, consensus, ref):
            new_labels = list()
            for label, consensus_letter, ref_letter in zip(labels, consensus, ref):
                if consensus_letter != ref_letter:
                    label = "[%s%s]" % (ref_letter, label)
                new_labels.append(label)
            return new_labels

        # make sequence labels
        consensus = self.consensus()
        if seqnum_labels:
            assert len(seqnum_labels) == self.max_length(), "Number of labels does not match sequence length"
        elif ref_id:
            refseq = self.get_seq_by_id(ref_id)
            if "seqnums" in refseq.letter_annotations:
                seqnum_labels = refseq.letter_annotations["seqnums"]
            # if reference sequence provided, highlight the consensus positions
            seqnum_labels = highlight_consensus(seqnum_labels, consensus, refseq)
        else:
            seqnum_labels = range(1, len(consensus) + 1)

        # set weblogo options
        opts = wl.LogoOptions(
            formatter=wl.formatters[logo_format],
            stacks_per_line=cols,
            color_scheme=color_scheme(),
            logo_font=logo_font,
            resolution=resolution,
            **kwargs,
        )
        opts.rotate_numbers = True
        opts.annotate = seqnum_labels

        # remove gap and stop characters so that they are not included in weblogo
        ignore = [gap_letter, stop_letter]
        alphabet = wl.seq.Alphabet("".join(s for s in self.alphabet if s not in ignore))
        counts = [count for s, count in self.as_counts_by_alphabet() if s not in ignore]

        # count position-specific frequencies
        data = wl.LogoData.from_counts(alphabet, np.array(counts).T)

        # return logo as png image
        logo = opts.formatter(data, wl.LogoFormat(data, opts))
        return Image.open(io.BytesIO(logo))

    def align(self, preserve_order: bool = True, *args, **kwargs):
        """Returns a Series of aligned SeqLikes from the specified column (by
        default, 'seqs').

        :param preserve_order: Whether or not to preserve the order of sequences.
        :param *args: Not used.
        :param **kwargs: Passed into the `seqlike.alignment_utils.align` function.
        :returns: A pandas Series of aligned sequences.
        """
        col_seq_type = self._obj.iloc[0]._type
        alignment = align(self._obj, seq_type=col_seq_type, preserve_order=preserve_order, **kwargs)
        return pd.Series([SeqLike(x, col_seq_type) for x in alignment], self._obj.index)

    def as_alignment(self, alphabet: Optional[str] = None) -> MultipleSeqAlignment:
        """Return a `Bio.Align.MultipleSeqAlignment` of the specified columns'
        SeqRecords.

        Because MSAs must be the same length, we pad the ends if
        needed.  This allows for plotting of SeqLikes of variable
        lengths.

        :param alphabet: The SeqLike alphabet to use.
        :returns: A `Bio.Align.MultipleSeqAlignment` object.
        """
        seq_lens = self._obj.apply(len).tolist()
        if len(set(seq_lens)) != 1:
            max_len = max(seq_lens)
            seqs = [x.pad_to(max_len) for x in self._obj]
        else:
            seqs = self._obj
        return MultipleSeqAlignment([x.to_seqrecord(alphabet=alphabet) for x in seqs])

    def as_counts(self, pad=True, dtype=float, encoder=None) -> np.ndarray:
        """Return a 2D numpy array of letter counts from the sequences or alignment.,

        Here the sequence position indices (alignment columns)
        are the columns of the array,
        and the rows correspond to the letters.

        :param pad: Whether or not to pad sequence.
        :param dtype: numpy dtype
        :param encoder: The one-hot encoder to use.
        :returns: A NumPy array.
        """
        if encoder is None:
            encoder = onehot_encoder_from_alphabet(self.alphabet)
        return self.to_onehot(pad=pad, dtype=dtype, encoder=encoder).sum(axis=0)

    def as_counts_df(self, pad=True, dtype=float, encoder=None):
        """Return DataFrame of letter counts from the sequences (or sequence alignment).

        The sequence position indices (alignment columns) are the columns of the array,
        and the rows correspond to the letters in the alphabet.

        :param pad: Whether or not to pad sequence.
        :param dtype: numpy dtype
        :param encoder: The one-hot encoder to use.
        :returns: A pandas DataFrame
        """
        return pd.DataFrame(
            self.as_counts(pad=pad, dtype=dtype, encoder=encoder),
            columns=list(self.alphabet),
        ).T

    def as_counts_by_alphabet(self, pad=True, dtype=float, encoder=None):
        """
        Return generator of (alphabet letter, letter counts) tuples.

        This is done for each letter in alphabet,
        where the letter counts are indexed by column of the sequence alignment.

        :param pad: Whether or not to pad sequence.
        :param dtype: numpy dtype
        :param encoder: The one-hot encoder to use.
        :returns: A generator of (alphabet letter, count) tuples
        """
        return zip(self.alphabet, self.as_counts(pad=pad, dtype=dtype, encoder=encoder).T)

    @property
    def alphabet(self) -> str:
        """Return the alphabet string of the all of the sequences.

        :returns: An alphabet string.
        """
        alphabets = self._obj.apply(lambda s: "".join(s.alphabet))
        assert len(alphabets.unique()) == 1
        return alphabets.iloc[0]

    def _extend_ambiguous_counts(self) -> np.ndarray:
        """Distribute ambiguous letter counts among represented unambiguous letters.
        :returns: a numpy 2d array of letter counts by sequence position after expanding
            the ambiguous letters.
        """
        if self._type == "AA":
            ambiguous_values = extended_protein_values
        else:
            ambiguous_values = ambiguous_nt_values
        # counts of letter identity by sequence position
        counts = self.as_counts()
        new_counts = np.zeros(counts.shape)
        alphabet = self.alphabet
        # expand each letter in the ambiguous alphabet
        for j, letter in enumerate(alphabet):
            if letter in ambiguous_values:
                for unambiguous_letter in ambiguous_values[letter]:
                    i = alphabet.index(unambiguous_letter)
                    # counts matrix has shape (len(seq), len(alphabet))
                    new_counts[:, i] += counts[:, j]
        return new_counts

    def consensus(self, ignore_gap=True) -> SeqLike:
        """Return the consensus sequence as a SeqLike.

        Ambiguous letter counts are distributed among represented unambiguous letters.

        :param ignore_gap: Whether to ignore gaps or not. Defaults to True.
        :returns: The consensus sequence as a SeqLike.
        """
        if self._type == "AA":
            alphabet = STANDARD_AA
        else:
            alphabet = STANDARD_NT
        counts = self._extend_ambiguous_counts()
        if ignore_gap:
            # zero the gap letter counts so that gaps do not show up in consensus
            counts[:, alphabet.index(gap_letter)] = 0
        sequence = ""
        for i in range(len(counts)):
            j = np.argmax(counts[i, :])
            if counts[i, j] > 0:
                sequence_letter = alphabet[j]
            else:
                # if and only if no consensus (zero non-gap counts) at this position, use gap letter
                sequence_letter = gap_letter
            sequence += sequence_letter
        return SeqLike(sequence, self._type)

    def degenerate(self) -> SeqLike:
        """Return the ambiguous sequence representation of the sequences as a SeqLike

        Following the rules adapted from
        D. R. Cavener: "Comparison of the consensus sequence flanking
        translational start sites in Drosophila and vertebrates."
        Nucleic Acids Research 15(4): 1353-1361. (1987).
        The same rules are used by TRANSFAC.

        :sa: http://biopython.org/DIST/docs/api/Bio.motifs.matrix-pysrc.html

        :returns: The degnerate SeqLike.
        """
        if self._type == "AA":
            alphabet = STANDARD_AA
            ambiguous_values = extended_protein_values
            generic_value = generic_protein_letter
        else:
            alphabet = STANDARD_NT
            ambiguous_values = ambiguous_nt_values
            generic_value = generic_nt_letter
        # ambiguous letter indexed by sorted representative unambiguous letters
        reverse_ambiguous_values = dict((v, k) for k, v in ambiguous_values.items())
        counts = self._extend_ambiguous_counts()
        sequence = ""
        for i in range(len(counts)):
            key = ""
            for j, letter in enumerate(alphabet):
                if counts[i, j] > 0 and letter is not gap_letter:
                    key += letter
            # if and only if no consensus (counts) at this position, use gap letter
            if len(key) == 0:
                sequence += gap_letter
            else:
                if len(key) > 1 and "TU" not in key:
                    if "T" in key:
                        key = key.replace("T", "TU")
                    elif "U" in key:
                        key = key.replace("U", "TU")
                try:
                    ambiguous_letter = reverse_ambiguous_values[key]
                except KeyError:
                    ambiguous_letter = generic_value
                sequence += ambiguous_letter
        return SeqLike(sequence, self._type)

    def __getitem__(self, index):
        """Slice sequences by dataframe row and sequence column.

        :param index: A tuple specifying the dataframe row and sequence column.
        :returns: The letter of a given sequence.
        """
        assert isinstance(index, tuple) and len(index), "Index is row and column"
        assert isinstance(index[0], (slice, int)) and isinstance(index[1], (slice, int, list))
        return self._obj.iloc[index[0]].apply(lambda seq: seq[index[1]])

    def max_length(self):
        return int(self._obj.apply(len).max())

    def get_seq_by_id(self, seq_id):
        """Get a sequence record by id.

        :param seq_id: The ID to search for.
        :returns: The first SeqLike object that has that particular ID.
        """
        seqrow = self._obj[self._obj.apply(lambda x: x.id == seq_id)].drop_duplicates()
        assert len(seqrow) == 1
        return seqrow.iloc[0]

    def slice_to_ref(self, ref_id, list_of_seqnums=None):
        """
        Slice alignment sequences by columns corresponding to seqnums.

        :param ref_id: The `id` of the sequence to use
            as a reference for position numbers.
        :param list_of_seqnums: A list of integers corresponding to the positions.
            This is optional; if not provided we defer to `seqnums` field
            in the `letter_annotations` of the SeqLike object.
        :returns: A pandas Series with the sliced SeqLikes.
        """
        # find the reference sequence
        refseq = self.get_seq_by_id(ref_id)
        # find the column indices corresponding to the seqnums
        if list_of_seqnums is None:
            list_of_seqnums = [
                seqnum for seqnum in refseq._seqrecord.letter_annotations["seqnums"] if seqnum is not None
            ]
        indices = refseq.seq_num_to_idx(list_of_seqnums)
        # slice the alignment at column indices
        return self.__getitem__((slice(None, None, None), indices))

    def nt(self, codon_map=None) -> pd.Series:
        """
        Return a Pandas Series of the NT form of the column of seqlikes.

        :returns: A pandas Series
        """
        return self._obj.apply(lambda x: x.nt(codon_map=codon_map))

    def aa(self) -> pd.Series:
        """
        Return a Pandas Series of the AA form of the column of seqlikes.

        :returns: A pandas Series
        """
        return self._obj.apply(lambda x: x.aa())

    def to_onehot(self, pad=True, dtype=float, encoder=None) -> np.ndarray:
        """Return a 3d Numpy array of the specified column in onehot encoding.

        The dimensions will be num_seqs x length x num_bases (5 for NT, 28 for AA)

        We pad if needed, because the numpy arrays must be the same size.

        :param pad: Whether or not to pad characters. Defaults to True.
        :param dtype: The dtype of the resulting numpy array.
        :param encoder: The sklearn-compatible encoder object.
            Defaults to None.
        :returns: A one-hot-encoded array.
        """
        if pad:
            max_len = self.max_length()
            return np.stack(
                self._obj.apply(lambda x: x.pad_to(max_len).to_onehot(dtype, encoder)).values,
                axis=0,
            )
        else:
            return np.stack(self._obj.apply(lambda x: x.to_onehot(dtype, encoder)).values, axis=0)

    def to_index(self, pad: bool = True, dtype: type = int, encoder=None) -> np.ndarray:
        """Return a 2d Numpy array of the specified column in index encoding.

        The dimensions will be num_seqs x length, and the values will
        range from 0 to num_bases-1 (4 for NT, 27 for AAs).

        We pad if needed, because the numpy arrays must be the same size.

        :param pad: Whether or not to pad characters. Defaults to True.
        :param dtype: The dtype of the resulting numpy array.
        :param encoder: The sklearn-compatible encoder object.
            Defaults to None.
        :returns: An index-encoded array.
        """
        if pad:
            max_len = self.max_length()
            return np.stack(
                self._obj.apply(lambda x: x.pad_to(max_len).to_index(dtype, encoder)).values,
                axis=0,
            )
        else:
            return np.stack(self._obj.apply(lambda x: x.to_index(dtype, encoder)).values, axis=0)

    def back_translate(self, codon_map=None):
        """Back-translate the collection of SeqLikes.

        :param codon_map: A SeqLike codon map to use.
        :returns: a Pandas Series of the specified column
            with back translated the AAs.
            Use the specified codon_map if given,
            or the codon_map in each SeqLike if not.
        """
        return self._obj.apply(lambda x: x.back_translate(codon_map=codon_map))

    def ungap(self):
        """Return ungapped seqlikes.

        Note that this may mean that this may disrupt any NT/AA correspondence.

        :returns: A Pandas series of the specified column with all gaps removed.
        """
        return self._obj.apply(lambda x: x.ungap())

    @property
    def _type(self):
        """Return a string that is the type of all SeqLikes in the 'seqs'
        column.

        Because this is a property, not a method, we can't pass in an optional column.

        :returns: The `_type` property of the _first_ seqlike object
            in the `seqs` column.
        """
        return self._obj.iloc[0]._type.upper()

    def __repr__(self):
        return f"{self._obj.__repr__()}"
