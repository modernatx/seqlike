from copy import deepcopy
from typing import List

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

from .alphabets import gap_letter, stop_letter
from .alignment_commands import mafft_alignment
from .SeqLike import SeqLikeType, SeqLike


def copy_annotations_from_unaligned(aligned_seqrec: SeqRecord, unaligned_seqrec: SeqRecord):
    """Copy letter annotations (and annotations) from unaligned to (aligned) seqrec.
    Assumes that seqrec is a gapped version of unaligned_seqrec

    :param aligned_seqrec: a SeqRecord from alignment
    :param unaligned_seqrec: unaligned SeqRecord corresponding to seqrec;
        contains reference letter_annotations
    :returns: a SeqRecord with letter_annotations that have been aligned
    """
    # NCBI Blast id includes description, whereas alignment does not
    assert aligned_seqrec.id in unaligned_seqrec.id, f"{aligned_seqrec.id} <> {unaligned_seqrec.id}"
    # copy annotations from previous
    newrec = deepcopy(aligned_seqrec)
    newrec.annotations = unaligned_seqrec.annotations
    # clear any letter annotations added during deepcopy
    newrec.letter_annotations = dict()
    # original sequence and letter annotations
    seq = unaligned_seqrec.seq
    letter_annotations = unaligned_seqrec.letter_annotations
    # index to track position in original sequence
    i = 0
    for j, letter in enumerate(aligned_seqrec.seq):
        if letter in [gap_letter, stop_letter]:
            for key, values in letter_annotations.items():
                # convert strings into lists of characters,
                # then combine into string at end of loop
                if key == "seqnums":
                    letter_annotation = None
                elif all(isinstance(value, str) for value in values):
                    letter_annotation = gap_letter
                else:
                    letter_annotation = None
                newrec.letter_annotations.setdefault(key, list()).append(letter_annotation)
        else:
            while seq[i] in [gap_letter, stop_letter]:
                i += 1
            assert letter == seq[i], f"letter {letter} at {j} <> seq {seq[i]} at {i}"
            for key in letter_annotations.keys():
                newrec.letter_annotations.setdefault(key, list()).append(letter_annotations[key][i])
            i += 1
    # convert list of chars into string
    for key, values in letter_annotations.items():
        if isinstance(values, str):
            newrec.letter_annotations[key] = "".join(newrec.letter_annotations[key])
    return newrec


def align_letter_annotations(unaligned, aligned, seq_type, original_ids=None):
    """Align letter_annotations from seqrecs to their alignment (aligned).
    Assumes that aligned is an alignment of unaligned.

    :param unaligned: iterator of SeqRecord
    :param aligned: iterator of SeqRecord, aligned version of unaligned
    :param seq_type: 'DNA', 'RNA', 'AA'
    :param original_ids: dict of original ID indexed by new ID; for proper ordering after alignment,
        seqids are changed to simple integer index; this argument helps restore the original IDs
    :returns: a new MultipleSeqAlignment with aligned letter_annotations
    """
    unaligned_dict = SeqIO.to_dict(unaligned)
    realigned = list()
    for seqrec in aligned:
        realigned.append(copy_annotations_from_unaligned(seqrec, unaligned_dict[seqrec.id]))
    if original_ids:
        # restore ID, name, and description; these will have changed upon readback from alignment output
        realigned = [SeqLike(s, seq_type=seq_type).to_seqrecord(**original_ids[s.id]) for s in realigned]
    return MultipleSeqAlignment(realigned)


def align(
    seqs: List[SeqLikeType],
    seq_type: str = "dna",
    aligner=mafft_alignment,
    preserve_order=True,
    **kwargs,
):
    """Align target sequences, preserving letter_annotations.

    This function maps original SeqRecord IDs to list of unique integer IDs for alignment,
    in case of common IDs, then replaces the original IDs when aligning letter_annotations

    :param seqs: an iterator of SeqRecord that will be aligned
    :param seq_type: 'DNA', 'RNA', 'AA'
    :param aligner: an alignment wrapper function from seqlike.alignment_commands
    :param preserve_order: if True, reorder aligned seqrecs to match input order.
    :param **kwargs: additional arguments for alignment command
    :returns: a MultipleSeqAlignment object with aligned sequences as SeqRecord
    """
    # map original IDs to unique IDs
    seqrecs = list()
    original_ids = dict()
    for i, seq in enumerate(seqs):
        seqlike = SeqLike(seq, seq_type=seq_type)
        # replace the ID
        seqrecs.append(seqlike.to_seqrecord(id=str(i)))
        # record original IDs, names, descriptions, as these will change upon readback from alignment output FASTA
        original_ids[str(i)] = {
            "id": seqlike.id,
            "name": seqlike.name,
            "description": seqlike.description,
        }

    # execute alignment with unique IDs
    aligned = aligner(seqrecs, preserve_order=preserve_order, **kwargs)

    # re-align the letter_annotations, original IDs/names/descriptions
    return align_letter_annotations(seqrecs, aligned, seq_type, original_ids)
