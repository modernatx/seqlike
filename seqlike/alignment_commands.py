# from __future__ import absolute_import, division, print_function, unicode_literals
from typing import List

import tempfile
from io import StringIO
import lazy_loader as lazy


from Bio import AlignIO, SeqIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import ClustalOmegaCommandline

from .AlignCommandline import MafftCommandline, MuscleCommandline
from .SeqLike import SeqLikeType, SeqLike

pd = lazy.load("pandas")


def pad_seq_records_for_alignment(seqs: List[SeqLikeType]):
    """Pad sequences so that lengths match for multiple sequence alignment.

    :param seqs: a list of SeqLikeType
    :returns: a MultipleSeqAlignment object
    """
    df = pd.DataFrame({"seqs": [SeqLike(seq, seq_type="aa") for seq in seqs]})
    return df.seqs.seq.as_alignment()


def _generic_aligner_commandline_stdout(cline, **kwargs):
    """Execute aligner commandline that writes to stdout and return an alignment. Helper function.

    :param cline: a subprocess object from Bio.Align.Applications.AbstractCommandline
    :param **kwargs: additional arguments for alignment command
    :returns: a MultipleSeqAlignment object
    """
    stdout, _ = cline()
    try:
        stdout = StringIO(stdout)
    except TypeError:
        stdout = StringIO(unicode(stdout, "utf-8"))
    return AlignIO.read(stdout, "fasta", **kwargs)


def _generic_aligner_commandline_file(cline, seqrecs, **kwargs):
    """Execute aligner commandline that requires file i/o and return an alignment. Helper function.

    :param cline: a subprocess object from Bio.Align.Applications.AbstractCommandline
    :param seqrecs: a list of SeqRecord that will be aligned
    :param **kwargs: additional arguments for alignment command
    :returns: a MultipleSeqAlignment object
    """
    assert len(seqrecs) > 1, "Need more than 1 sequence for alignment."
    # build alignment object 'unaligned'; pad seqrecs to be equal length
    unaligned = pad_seq_records_for_alignment(seqrecs)
    # execute alignment
    with tempfile.NamedTemporaryFile(delete=False, mode="w") as tempf:
        AlignIO.write(unaligned, tempf, "fasta")
        tempf.flush()
        return cline(tempf, **kwargs)


def _generic_alignment(cline, seqrecs, preserve_order=True, **kwargs):
    """Align sequences using command line stored as cline. Helper function.

    :param cline: a subprocess object from Bio.Align.Applications.AbstractCommandline
    :param seqrecs: an iterator of SeqRecord that will be aligned
    :param preserve_order: if True, reorder aligned seqrecs to match input order.
    :param **kwargs: additional arguments for alignment command
    :returns: a MultipleSeqAlignment object with aligned sequences
    """
    # convert iterator to list, so that we can extract keys and still run the alignment
    unaligned = list(seqrecs)
    # if alignment sequences from NCBI Blast, id will include spaces
    keys = [seqrec.id.split()[0] for seqrec in unaligned]
    # execute alignment
    aligned = _generic_aligner_commandline_file(cline, unaligned, **kwargs)
    if preserve_order:
        aligned = SeqIO.to_dict(aligned)
        aligned = MultipleSeqAlignment(aligned[key] for key in keys)
    # make all alignment uppercase
    return MultipleSeqAlignment([seqrec.upper() for seqrec in aligned])


def mafft_alignment(seqrecs, preserve_order=True, **kwargs):
    """Align sequences using MAFFT.

    :param seqrecs: a list or dict of SeqRecord that will be aligned to ref
    :param preserve_order: if True, reorder aligned seqrecs to match input order.
    :param **kwargs: additional arguments for alignment command
    :returns: a MultipleSeqAlignment object with aligned sequences

    :sa: https://mafft.cbrc.jp/alignment/software/
    """

    def commandline(file_obj, **kwargs):
        cline = MafftCommandline(input=file_obj.name, **kwargs)
        return _generic_aligner_commandline_stdout(cline)

    # MAFFT does not reorder alignment by default (reorder=False), but don't overwrite 'reorder' if set
    if "reorder" not in kwargs:
        kwargs["reorder"] = not preserve_order
    return _generic_alignment(commandline, seqrecs, preserve_order=preserve_order, **kwargs)


def muscle_alignment(seqrecs, preserve_order=True, **kwargs):
    """Align sequences using Muscle 3.8.

    Notes:

    1. Muscle's latest version is version 5.1. However, the interface is different from version 3.8
       BioPython's MuscleCommandline is only compatible with version 3.8.
       See this comment for more information: https://github.com/modernatx/seqlike/pull/60#issue-1257542946
    2. The preserve_order parameter (preserves original sequence order,
       as aligner may try to group sequences by similarity) may still be buggy.
       Also see this comment for more information: https://github.com/modernatx/seqlike/pull/60#issue-1257542946

    :param seqrecs: a list or dict of SeqRecord that will be aligned to ref
    :param preserve_order: if True, reorder aligned seqrecs to match input order.
    :param **kwargs: additional arguments for alignment command
    :returns: a MultipleSeqAlignment object with aligned sequences
    """

    def commandline(file_obj, **kwargs):
        cline = MuscleCommandline(input=file_obj.name, **kwargs)
        return _generic_aligner_commandline_stdout(cline)

    # Muscle reorders alignment by default, but don't overwrite 'group' if already set
    if "group" not in kwargs:
        kwargs["group"] = not preserve_order
    return _generic_alignment(commandline, seqrecs, preserve_order=preserve_order, **kwargs)


def clustal_omega_alignment(seqrecs, preserve_order=True, **kwargs):
    """Align sequences using Clustal Omega

    :param seqrecs: a list or dict of SeqRecord that will be aligned to ref
    :param preserve_order: if True, reorder aligned seqrecs to match input order.
    :param **kwargs: additional arguments for alignment command
    :returns: a MultipleSeqAlignment object with aligned sequences
    """
    if preserve_order:
        outputorder = "input-order"
    else:
        outputorder = "tree-order"

    def commandline(file_obj, **kwargs):
        cline = ClustalOmegaCommandline("clustalo", infile=file_obj.name, outputorder=outputorder, **kwargs)
        return _generic_aligner_commandline_stdout(cline)

    return _generic_alignment(commandline, seqrecs, **kwargs)


def clustal_omega_distance_matrix(seqrecs, **kwargs):
    """Generate a distance matrix using Clustal Omega

    :param seqrecs: a list or dict of SeqRecord that will be aligned to ref
    :param **kwargs: additional arguments for command line alignment
    :returns: the pairwise distance matrix
    """

    def commandline(ft, **kwargs):
        with tempfile.NamedTemporaryFile(delete=False, mode="w") as ft_out:
            cline = ClustalOmegaCommandline(
                "clustalo",
                infile=ft.name,
                force=True,
                distmat_out=ft_out.name,
                distmat_full=True,
                distmat_full_iter=True,
            )
        stdout, stderr = cline()
        df = pd.read_csv(ft_out.name, delim_whitespace=True, skiprows=1, header=None, index_col=0)
        df.index.name = "seqid"
        return df

    return _generic_aligner_commandline_file(commandline, seqrecs, **kwargs)


def clustal_omega_alignment_tree(seqrecs, **kwargs):
    """Generate phylogenetic tree using Clustal Omega and scikit-bio Neighbor Joining

    This function computes a distance matrix using Clustal Omega, which skbio.tree.nj
    uses to generate a newick file. Bio.Phylo can read this newick file.

    Note: this function requires scikit-bio.

    :param seqrecs: a list or dict of SeqRecord that will be aligned to ref
    :param **kwargs: additional arguments for alignment command
    :returns: a Bio.Phylo phylogenetic tree object

    :sa: https://biopython.org/wiki/Phylo
    :sa: http://scikit-bio.org/docs/0.2.1/generated/skbio.tree.nj.html
    """
    import skbio

    def skbio2phylo(treenode, format="newick"):
        """Convert skbio.tree.TreeNode object to Bio.Phylo.Newick.Tree object

        :param treenode: an skbio.tree.TreeNode object
        :param format: kind of tree, AKA New Hampshire Format
        :returns: an equivalent Bio.Phylo.Newick.Tree object

        :sa: https://biopython.org/docs/1.74/api/Bio.Phylo.Newick.html
        """
        with tempfile.NamedTemporaryFile(delete=True, mode="w") as tempf:
            treenode.write(tempf.name, format)
            tempf.flush()
            return Phylo.read(tempf.name, format)

    distance_matrix = clustal_omega_distance_matrix(seqrecs, **kwargs)
    ids = [s.id for s in seqrecs]
    skbio_tree = skbio.tree.nj(skbio.DistanceMatrix(distance_matrix, ids))
    return skbio2phylo(skbio_tree)


def clustalw_alignment_tree(seqrecs, **kwargs):
    """Generate phylogenetic tree using ClustalW. Note that ClustalW is an older
    generation of Clustal aligner compared to Clustal Omega. It is considered to
    be slower and less robust for aligning large sequence sets, but it is included here
    because it does not require scikit-bio.

    :param seqrecs: a list or dict of SeqRecord that will be aligned to ref
    :param **kwargs: additional arguments for alignment command
    :returns: the phylogenetic tree instead of the alignment object
    """

    def commandline(ft, **kwargs):
        with tempfile.NamedTemporaryFile(delete=False, mode="w") as ft_out:
            cline = ClustalwCommandline(infile=ft.name, output="fasta", newtree=ft_out.name)
            stdout, stderr = cline()
            return Phylo.read(ft_out.name, "newick")

    return _generic_alignment(commandline, seqrecs, preserve_order=False, **kwargs)
