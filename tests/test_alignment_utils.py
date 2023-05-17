from Bio.SeqRecord import SeqRecord
from seqlike.alignment_utils import copy_annotations_from_unaligned


def test_copy_annotations_from_unaligned():
    aligned_seqstr = "MFA-TS*"
    aligned_seqnums = ["1", "2", "3", None, "4", "5", "6"]
    unaligned_seqstr = "MFATS*"
    unaligned_seqnums = ["1", "2", "3", "4", "5", "6"]

    aligned_seqrec = SeqRecord(aligned_seqstr, id="aligned")
    # aligned_seqrec.letter_annotations = {"seqnums": aligned_seqnums}

    unaligned_seqrec = SeqRecord(unaligned_seqstr, id="unaligned")
    unaligned_seqrec.letter_annotations = {"seqnums": unaligned_seqnums}

    new_aligned_seqrec = copy_annotations_from_unaligned(aligned_seqrec, unaligned_seqrec)
    assert str(new_aligned_seqrec.seq) == str(aligned_seqrec.seq) == aligned_seqstr
    assert new_aligned_seqrec.letter_annotations["seqnums"] == aligned_seqnums
