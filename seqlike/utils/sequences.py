from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..alphabets import gap_letter


def add_seqnums_to_letter_annotations(seqrec, include_gaps=True):
    # If seqnums is already present, don't override it.
    if "seqnums" in seqrec.letter_annotations:
        return seqrec
    if include_gaps:
        seqnums = [str(i + 1) for i in range(len(seqrec))]
    else:
        seqnums = list()
        for aa in seqrec.seq:
            if aa == gap_letter:
                seqnums.append(None)
            else:
                seqnums.append(str(len(seqnums) + 1))
    seqrec.letter_annotations["seqnums"] = seqnums
    return seqrec


def ungap(seqlike, gap=None):
    """Return a SeqRecord without the gap character(s) in the sequence.

    The gap character can be specified in two ways -
    either as an explicit argument,
    or via the sequence's alphabet.

    For example:

    ```python
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna

    my_dna = SeqRecord(Seq("-ATA--TGAAAT-TTGAAAA-", generic_dna), id="X")
    my_dna.seq
    # Output: Seq('-ATA--TGAAAT-TTGAAAA-', DNAAlphabet())

    my_dna.ungap("-").seq
    # Output: Seq('ATATGAAATTTGAAAA', DNAAlphabet())
    ```

    If the gap character is not given as an argument,
    it will be taken from the sequence's alphabet (if defined).
    For more details, see the Seq object's ungap method.

    Any per-letter-annotation is sliced
    to match how the sequence gets sliced to remove the gaps.
    Other annotations is retained as is.
    SeqFeature locations are adjusted to use the new coordinates:

    ```python
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    my_dna.features.append(SeqFeature(FeatureLocation(0,4)))
    my_dna.features.append(SeqFeature(FeatureLocation(0,5)))
    my_dna.features.append(SeqFeature(FeatureLocation(1,4)))
    my_dna.features.append(SeqFeature(FeatureLocation(1,12)))
    my_dna.features.append(SeqFeature(FeatureLocation(4,13)))
    my_dna.features.append(SeqFeature(FeatureLocation(5,13)))
    my_dna.features.append(SeqFeature(FeatureLocation(6,13)))
    my_dna.features.append(SeqFeature(FeatureLocation(12,20)))
    my_dna.features.append(SeqFeature(FeatureLocation(12,21)))
    for f in my_dna.features:
        print f.location, f.extract(my_dna.seq)
    ```

    ```
    # The output is below
    [0:4] -ATA
    [0:5] -ATA-
    [1:4] ATA
    [1:12] ATA--TGAAAT
    [4:13] --TGAAAT-
    [5:13] -TGAAAT-
    [6:13] TGAAAT-
    [12:20] -TTGAAAA
    [12:21] -TTGAAAA-
    ```

    Notice most of these examples deliberately have
    the features ending on a gap character.
    The start and end positions are adjusted to ensure
    the feature describes the equivalent ungapped sequence:

    ```python
    ungapped = my_dna.ungap("-")
    for f in ungapped.features:
        print(f.location, f.extract(ungapped.seq))
    ```

    ```
    # Output
    [0:3] ATA
    [0:3] ATA
    [0:3] ATA
    [0:9] ATATGAAAT
    [3:9] TGAAAT
    [3:9] TGAAAT
    [3:9] TGAAAT
    [9:16] TTGAAAA
    [9:16] TTGAAAA
    ```

    For example with per-letter-annotation, we'll use Bio.SeqIO to load an
    Ace assembly which includes quality scores but the sequence will be
    padded with any gap characters (for which there is no quality score
    available). You may want to get the ungapped sequence with its quality
    scores (e.g. to output as FASTQ):

    ```python
    from Bio import SeqIO
    record = SeqIO.read("Ace/consed_sample.ace", "ace")
    print(len(record))
    # 1475
    print(len(record) - record.seq.count("-"))
    # 1468
    print(record[860:880].format("fastq"))
    # @Contig1
    # CAGCAGAGAAGGGTTTGAAA
    # +
    # z{{{{yyyyyy{{{{{{{{{
    # <BLANKLINE>
    ```
    In the above example we've selected a subsection of the record to show
    in FASTQ format. Now lets remove the gaps:

    ```python
    ungapped = record.ungap()
    print(len(ungapped))
    # 1468
    ```

    Notice below how the coordinates for the region [860:880] have shifted
    by two since there are two gaps before it in the original record:

    ```python
    record.seq[0:860].count("-")
    # 2
    record[860:880].format("fastq") == ungapped[858:878].format("fastq")
    # True
    ```

    So, using this method we can take the gapped consensus records from any
    ACE file and save them as ungapped records in FASTQ, FASTA, QUAL, etc:

    ```python
    records = (rec.ungap() for rec in SeqIO.parse(in_file, "ace"))
    count = SeqIO.write(records, out_file, "fastq")
    ```

    Ungap method by Peter Cock
    From https://github.com/peterjc/biopython/blob/ace-reads/Bio/SeqRecord.py

    Based on discussion: http://www.biopython.org/pipermail/biopython-dev/2010-June/007878.html

    ## Historical Changelog

    31-Mar-2016: copied here by adousis
    26-Apr-2016: split 'ungap' into 'slice' helper function and 'ungap'
    24-Jun-2019: copied out of BioExtensions by agiessel

    :param seqlike: The seqlike object to ungap.
    :param gap: The gap character to remove from the sequence.
        If not specified, then the SeqLike alphabet's character is used.
    :returns: Ungapped seqlike.
    """
    if not gap:
        try:
            gap = seqlike.seq.alphabet.gap_char
        except AttributeError:
            gap = gap_letter
    new_seq = seqlike.seq.ungap(gap)

    if str(new_seq) == str(seqlike.seq):
        # Unchanged, not even the alphabet - don't need to alter annotation
        return seqlike
    if not gap:
        gap = seqlike.seq.alphabet.gap_char

    if not (seqlike.features or seqlike.letter_annotations):
        return SeqRecord(
            new_seq,
            id=seqlike.id,
            name=seqlike.name,
            description=seqlike.description,
            dbxrefs=seqlike.dbxrefs[:],
            annotations=seqlike.annotations.copy(),
        )

    slices = []
    new_index = -1
    if seqlike.seq[0] == gap:
        start = None
        in_gap = True
    else:
        start = 0
        in_gap = False
    mapping = []
    for old_index, letter in enumerate(seqlike):
        if letter == gap:
            mapping.append(None)
            if in_gap:
                pass
            else:
                in_gap = True
                if start is not None:
                    slices.append((start, old_index))
                    start = None
        else:
            new_index += 1
            assert letter == new_seq[new_index]
            mapping.append(new_index)
            if in_gap:
                in_gap = False
                assert start is None
                start = old_index
    if not in_gap:
        slices.append((start, len(seqlike)))
    assert len(mapping) == len(seqlike)
    mapping.append(len(new_seq))
    answer = slice_seqrec(seqlike, slices, ref_seq=new_seq)
    # Apply the mapping to the feature coordinates
    for f in seqlike.features:
        answer.features.append(f._map(mapping))
    return answer


def slice_seqrec(seqlike, slices, ref_seq=None):
    def slice_str(s, slices):
        return "".join(s[start:end] for start, end in slices)

    new_seq = Seq(slice_str(str(seqlike.seq), slices))
    if ref_seq:
        if str(ref_seq) != new_seq:
            msg = "%s\n%s\n" % (repr(slices), seqlike.seq)
            for s, e in slices:
                msg += " " * s + seqlike.seq[s:e] + "\n"
            assert False, msg
        new_seq = ref_seq

    answer = SeqRecord(
        new_seq,
        id=seqlike.id,
        name=seqlike.name,
        description=seqlike.description,
        dbxrefs=seqlike.dbxrefs[:],
        annotations=seqlike.annotations.copy(),
    )
    # Apply the slices to the per letter annotation
    for key, value in seqlike.letter_annotations.items():
        s, e = slices[0]
        new = value[s:e]
        for s, e in slices[1:]:
            new += value[s:e]
        answer.letter_annotations[key] = new
    return answer
