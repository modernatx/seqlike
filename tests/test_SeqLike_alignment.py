from types import *
import pandas as pd
from seqlike import SeqLike
from seqlike.alignment_utils import align
from seqlike.alignment_commands import mafft_alignment, muscle_alignment, pad_seq_records_for_alignment
import pytest


def get_aa_seqrecs():
    seqs = [
        ("P0", "KTQSVLYNNKYDYVWGSASTGHYRYNYFDY"),
        ("P1", "KTQSVIYNNKYDYVWGSYRYNYFDY"),
        ("P2", "KTQSVLYYDYVYGSASTGHYRYNHFDY"),
        ("P3", "SVLYNQKYDYVWGSASTGHYRYNFDY"),
        ("P4", "KTQSVNNKYDYVWGSASTGHYRFDY"),
    ]
    return [SeqLike(seq, seq_type="aa", id=seq_id).to_seqrecord() for seq_id, seq in seqs]


def get_aligned_aa_seqrecs():
    seqs = [
        ("P0", "KTQSVLYNNKYDYVWGSASTGHYRYNYFDY"),
        ("P1", "KTQSVIYNNKYDYVWGS-----YRYNYFDY"),
        ("P2", "KTQSVLY---YDYVYGSASTGHYRYNHFDY"),
        ("P3", "---SVLYNQKYDYVWGSASTGHYRYN-FDY"),
        ("P4", "KTQSV--NNKYDYVWGSASTGHYR---FDY"),
    ]
    return [SeqLike(seq, seq_type="aa", id=seq_id).to_seqrecord() for seq_id, seq in seqs]


def get_aligned_aa_seqrecs_duplicate_ids():
    seqs = [
        ("P0", "KTQSVLYNNKYDYVWGSASTGHYRYNYFDY"),
        ("P1", "KTQSVIYNNKYDYVWGS-----YRYNYFDY"),
        ("P2", "KTQSVLY---YDYVYGSASTGHYRYNHFDY"),
        # duplicate seq and ID
        ("P2", "KTQSVLY---YDYVYGSASTGHYRYNHFDY"),
        ("P3", "---SVLYNQKYDYVWGSASTGHYRYN-FDY"),
        # dubplicate ID, not duplicate seq
        ("P3", "KT-SVLYNQKYDYVWGSASTGHYRYN-FDY"),
        ("P4", "KTQSV--NNKYDYVWGSASTGHYR---FDY"),
    ]
    return [SeqLike(seq, seq_type="aa", id=seq_id).to_seqrecord() for seq_id, seq in seqs]


def get_nt_seqrecs():
    seqs = [
        ("A00", "ATGAGAGATTCACCATAA", [50] * 18),
        ("A01", "TTCACCATAA", [50] * 10),
        ("A02", "ATGAGAGATTCACCATAA", [50] * 18),
        ("B01", "TTCACAATAA", [50] * 10),
        ("B02", "ATGGGAGATTCACCATAA", [50] * 18),
    ]
    seqrecs = list()
    for (seq_id, seq, qual) in seqs:
        seqrec = SeqLike(seq, seq_type="dna", id=seq_id).to_seqrecord()
        seqrec.letter_annotations["phred_quality"] = qual
        seqrecs.append(seqrec)
    return seqrecs


@pytest.mark.xfail(reason="May fail if MAFFT is not installed.")
def test_mafft_alignment():
    seqrecs = get_aa_seqrecs()
    # test mafft (but not really testing anything except that the command works...)
    aligned1 = mafft_alignment(seqrecs)
    aligned2 = mafft_alignment(seqrecs, dash=True)


@pytest.mark.xfail(reason="May fail if Muscle is not installed.")
def test_muscle_alignment():
    seqrecs = get_aa_seqrecs()
    # test muscle
    aligned1 = muscle_alignment(seqrecs)
    aligned2 = muscle_alignment(seqrecs, group=False)
    assert seqrecs[1].id == aligned2[1].id
    assert seqrecs[1].seq != aligned2[1].seq
    assert aligned1 != aligned2


def test_pad_seq_records_for_alignment():
    nt_seqrecs = get_nt_seqrecs()
    padded_nt_seqrecs = pad_seq_records_for_alignment(nt_seqrecs)
    assert len(padded_nt_seqrecs[0]) == len(nt_seqrecs[0])

    aa_seqrecs = get_aa_seqrecs()
    padded_aa_seqrecs = pad_seq_records_for_alignment(aa_seqrecs)
    for seqrec in padded_aa_seqrecs[1:]:
        assert len(padded_aa_seqrecs[0]) == len(seqrec)


def test_align_nucleotide():
    seqrecs = get_nt_seqrecs()
    aligned = align(seqrecs)

    # test letter_annotations
    assert aligned[0].letter_annotations["phred_quality"] == [50] * 18
    assert aligned[1].letter_annotations["phred_quality"] == [None] * 8 + [50] * 10
    assert aligned[2].letter_annotations["phred_quality"] == [50] * 18
    assert aligned[3].letter_annotations["phred_quality"] == [None] * 8 + [50] * 10
    assert aligned[4].letter_annotations["phred_quality"] == [50] * 18


def test_align_protein():
    seqrecs = get_aa_seqrecs()
    ref_aligned = get_aligned_aa_seqrecs()

    # test basic alignment
    aligned = align(seqrecs, seq_type="aa")
    for aligned_seqrec, ref_seqrec in zip(aligned, ref_aligned):
        assert aligned_seqrec.id == ref_seqrec.id
        assert str(aligned_seqrec.seq) == str(ref_seqrec.seq)

    # test letter_annotations
    assert aligned[0].letter_annotations["seqnums"] == [str(i) for i in range(1, 31)]
    assert aligned[1].letter_annotations["seqnums"] == [str(i) for i in range(1, 18)] + [None] * 5 + [
        str(i) for i in range(18, 26)
    ]
    assert aligned[2].letter_annotations["seqnums"] == [str(i) for i in range(1, 8)] + [None] * 3 + [
        str(i) for i in range(8, 28)
    ]

    # test slice by column index
    col_indices = [4, 14, 16, 18]
    df = pd.DataFrame({"seqs": [SeqLike(seq, "aa") for seq in aligned]})
    sub_aligned = df.seqs.seq[:, col_indices]
    assert isinstance(sub_aligned, pd.Series)
    assert len(sub_aligned) == 5
    assert "".join(sub_aligned[0]) == "VWSS"
    assert "".join(sub_aligned[1]) == "VWS-"
    assert "".join(sub_aligned[2]) == "VYSS"

    # slice by seqnum
    seqnums = ["5", "15", "17", "19"]
    sub_aligned = df.seqs.seq.slice_to_ref("P0", seqnums)
    assert isinstance(sub_aligned, pd.Series)
    assert len(sub_aligned) == 5
    assert "".join(sub_aligned[0]) == "VWSS"
    assert "".join(sub_aligned[1]) == "VWS-"
    assert "".join(sub_aligned[2]) == "VYSS"
