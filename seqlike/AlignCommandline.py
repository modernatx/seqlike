"""
This extends Bio.Align
"""
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import *  # NOQA

from Bio.Align import Applications
from Bio.Application import AbstractCommandline, _Switch, _Option, _Argument


class BowtieCommandline(AbstractCommandline):
    """
    Interface to Bowtie2 alignment.

    The following example:

    ```python
    cline = BowtieCommandline(
        local=True,
        threads=4,
        extends=20,
        reseeds=3,
        mismatches=1,
        length=15,
        interval='S,1,0.5',
        index='ref',
        read1='read1.fastq',
        read2='read2.fastq',
        sam='test.sam'
    )
    stdout,stderr = cline()
    ```

    is equivalent to the command line:

    ```bash
    bowtie2 --local -p 4 -D 20 -R 3 -N 1 -L 15 -i S,1,0.5 -x ref -1 read1.fastq -2 read2.fastq -S test.sam
    ```
    """

    def __init__(self, cmd="bowtie2", **kwargs):
        self.program_name = cmd
        self.parameters = [
            _Switch(["--local", "local"], "perform local read alignment"),
            _Switch(["-r", "raw"], "query input files are raw one-sequence-per-line"),
            ## performance
            _Option(["-p", "threads"], "number of alignment threads to launch", equate=False),
            ## effort
            _Option(
                ["-D", "extends"],
                "give up after this many failed extends in a row",
                equate=False,
            ),
            _Option(
                ["-R", "reseeds"],
                'max # "re-seed" reads with repetitive seeds',
                equate=False,
            ),
            ## alignment
            _Option(
                ["-N", "mismatches"],
                "max # mismatches in seed alignment",
                checker_function=lambda x: x in [0, 1],
                equate=False,
            ),
            _Option(
                ["-L", "length"],
                "length of seed substring",
                checker_function=lambda x: x > 3,
                equate=False,
            ),
            _Option(
                ["-i", "interval"],
                "function for interval between seed substrings w/r/t read len",
                equate=False,
            ),
            ## main arguments
            _Option(["-x", "index"], "index filename prefix", equate=False),
            _Option(["-1", "read1"], "files with #1 mates", filename=True, equate=False),
            _Option(["-2", "read2"], "files with #2 mates", filename=True, equate=False),
            _Option(["-S", "sam"], "file for SAM output", filename=True, equate=False),
        ]
        super(BowtieCommandline, self).__init__(cmd, **kwargs)


class BowtieBuildCommandline(AbstractCommandline):
    """
    Interface to Bowtie2-build indexing software
    """

    def __init__(self, cmd="bowtie2-build", **kwargs):
        self.program_name = cmd
        self.parameters = [
            _Argument(["infile"], "input filename", filename=True),
            _Argument(["bt2_base"], "bt2 output basename"),
            _Switch(["-q", "--quiet", "quiet"], "print only error messages"),
        ]
        super(BowtieBuildCommandline, self).__init__(cmd, **kwargs)


class MafftCommandline(AbstractCommandline):
    def __init__(self, cmd="mafft", **kwargs):
        self.parameters = [
            # MAFFT-DASH is a protocol that utilizes structural alignment information available
            # in the DASH database to improve the MAFFT sequence alignment.
            _Switch(
                ["--dash", "dash"],
                "Use structural alignment information from DASH database.",
            ),
        ]
        # append the original set of parameters from MafftCommandline
        self.parameters += Applications.MafftCommandline().parameters
        AbstractCommandline.__init__(self, cmd, **kwargs)


class MuscleCommandline(AbstractCommandline):
    """
    Extend Biopython's MuscleCommandline to include two additional options

    \sa Bio.Align.Applications.MuscleCommandline
    """

    def __init__(self, cmd="muscle", **kwargs):
        ## This is taken from Applications.MuscleCommandline with some amendments
        CLUSTERING_ALGORITHMS = ["upgma", "upgmb", "neighborjoining"]
        DISTANCE_MEASURES_ITER1 = [
            "kmer6_6",
            "kmer20_3",
            "kmer20_4",
            "kbit20_3",
            "kmer4_6",
        ]
        DISTANCE_MEASURES_ITER2 = DISTANCE_MEASURES_ITER1 + [
            "pctid_kimura",
            "pctid_log",
        ]
        OBJECTIVE_SCORES = ["sp", "ps", "dp", "xp", "spf", "spm"]
        TREE_ROOT_METHODS = ["pseudo", "midlongestspan", "minavgleafdist"]
        SEQUENCE_TYPES = ["protein", "nucleo", "auto"]
        WEIGHTING_SCHEMES = [
            "none",
            "clustalw",
            "henikoff",
            "henikoffpb",
            "gsc",
            "threeway",
        ]
        self.parameters = [
            # Can't use "in" as the final alias as this is a reserved word in python:
            _Option(["-in", "in", "input"], "Input filename", filename=True, equate=False),
            _Option(["-out", "out"], "Output filename", filename=True, equate=False),
            _Switch(["-diags", "diags"], "Find diagonals (faster for similar sequences)"),
            _Switch(["-profile", "profile"], "Perform a profile alignment"),
            _Option(
                ["-in1", "in1"],
                "First input filename for profile alignment",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-in2", "in2"],
                "Second input filename for a profile alignment",
                filename=True,
                equate=False,
            ),
            # anchorspacing   Integer              32                 Minimum spacing between
            _Option(
                ["-anchorspacing", "anchorspacing"],
                "Minimum spacing between anchor columns",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # center          Floating point       [1]                Center parameter.
            #                                                        Should be negative.
            _Option(
                ["-center", "center"],
                "Center parameter - should be negative",
                checker_function=lambda x: isinstance(x, float),
                equate=False,
            ),
            # cluster1        upgma                upgmb              Clustering method.
            _Option(
                ["-cluster1", "cluster1"],
                "Clustering method used in iteration 1",
                checker_function=lambda x: x in CLUSTERING_ALGORITHMS,
                equate=False,
            ),
            # cluster2        upgmb                                   cluster1 is used in
            #                neighborjoining                         iteration 1 and 2,
            #                                                        cluster2 in later
            #                                                        iterations.
            _Option(
                ["-cluster2", "cluster2"],
                "Clustering method used in iteration 2",
                checker_function=lambda x: x in CLUSTERING_ALGORITHMS,
                equate=False,
            ),
            # diaglength      Integer              24                 Minimum length of
            #                                                        diagonal.
            _Option(
                ["-diaglength", "diaglength"],
                "Minimum length of diagonal",
                checker_function=lambda x: isinstance(x, int),
                equate=True,
            ),
            # diagmargin      Integer              5                  Discard this many
            #                                                        positions at ends of
            #                                                        diagonal.
            _Option(
                ["-diagmargin", "diagmargin"],
                "Discard this many positions at ends of diagonal",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # distance1       kmer6_6              Kmer6_6 (amino) or Distance measure for
            #                kmer20_3             Kmer4_6 (nucleo)   iteration 1.
            #                kmer20_4
            #                kbit20_3
            #                kmer4_6
            _Option(
                ["-distance1", "distance1"],
                "Distance measure for iteration 1",
                checker_function=lambda x: x in DISTANCE_MEASURES_ITER1,
                equate=False,
            ),
            # distance2       kmer6_6              pctid_kimura       Distance measure for
            #                kmer20_3                                iterations 2, 3 ...
            #                kmer20_4
            #                kbit20_3
            #                pctid_kimura
            #                pctid_log
            _Option(
                ["-distance2", "distance2"],
                "Distance measure for iteration 2",
                checker_function=lambda x: x in DISTANCE_MEASURES_ITER2,
                equate=False,
            ),
            # gapopen         Floating point       [1]                The gap open score.
            #                                                        Must be negative.
            _Option(
                ["-gapopen", "gapopen"],
                "Gap open score - negative number",
                checker_function=lambda x: isinstance(x, float),
                equate=False,
            ),
            # hydro           Integer              5                  Window size for
            #                                                        determining whether a
            #                                                        region is hydrophobic.
            _Option(
                ["-hydro", "hydro"],
                "Window size for hydrophobic region",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # hydrofactor     Floating point       1.2                Multiplier for gap
            #                                                        open/close penalties in
            #                                                        hydrophobic regions.
            _Option(
                ["-hydrofactor", "hydrofactor"],
                "Multiplier for gap penalties in hydrophobic regions",
                checker_function=lambda x: isinstance(x, float),
                equate=False,
            ),
            # log             File name            None.              Log file name (delete
            #                                                        existing file).
            _Option(["-log", "log"], "Log file name", filename=True, equate=False),
            # loga            File name            None.              Log file name (append
            #                                                        to existing file).
            _Option(
                ["-loga", "loga"],
                "Log file name (append to existing file)",
                filename=True,
                equate=False,
            ),
            # maxdiagbreak    Integer              1                  Maximum distance
            #                                                        between two diagonals
            #                                                        that allows them to
            #                                                        merge into one
            #                                                        diagonal.
            _Option(
                ["-maxdiagbreak", "maxdiagbreak"],
                "Maximum distance between two diagonals that allows " "them to merge into one diagonal",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # maxhours        Floating point       None.              Maximum time to run in
            #                                                        hours. The actual time
            #                                                        may exceed the
            #                                                        requested limit by a
            #                                                        few minutes. Decimals
            #                                                        are allowed, so 1.5
            #                                                        means one hour and 30
            #                                                        minutes.
            _Option(
                ["-maxhours", "maxhours"],
                "Maximum time to run in hours",
                checker_function=lambda x: isinstance(x, float),
                equate=False,
            ),
            # maxiters        Integer 1, 2 ...     16                 Maximum number of
            #                                                        iterations.
            _Option(
                ["-maxiters", "maxiters"],
                "Maximum number of iterations",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # maxtrees        Integer              1                  Maximum number of new
            #                                                        trees to build in
            #                                                        iteration 2.
            _Option(
                ["-maxtrees", "maxtrees"],
                "Maximum number of trees to build in iteration 2",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # minbestcolscore Floating point       [1]                Minimum score a column
            #                                                        must have to be an
            #                                                        anchor.
            _Option(
                ["-minbestcolscore", "minbestcolscore"],
                "Minimum score a column must have to be an anchor",
                checker_function=lambda x: isinstance(x, float),
                equate=False,
            ),
            # minsmoothscore  Floating point       [1]                Minimum smoothed score
            #                                                        a column must have to
            #                                                        be an anchor.
            _Option(
                ["-minsmoothscore", "minsmoothscore"],
                "Minimum smoothed score a column must have to " "be an anchor",
                checker_function=lambda x: isinstance(x, float),
                equate=False,
            ),
            # objscore        sp                   spm                Objective score used by
            #                ps                                      tree dependent
            #                dp                                      refinement.
            #                xp                                      sp=sum-of-pairs score.
            #                spf                                     spf=sum-of-pairs score
            #                spm                                     (dimer approximation)
            #                                                        spm=sp for < 100 seqs,
            #                                                        otherwise spf
            #                                                        dp=dynamic programming
            #                                                        score.
            #                                                        ps=average profile-
            #                                                        sequence score.
            #                                                        xp=cross profile score.
            _Option(
                ["-objscore", "objscore"],
                "Objective score used by tree dependent refinement",
                checker_function=lambda x: x in OBJECTIVE_SCORES,
                equate=False,
            ),
            # root1           pseudo               pseudo             Method used to root
            _Option(
                ["-root1", "root1"],
                "Method used to root tree in iteration 1",
                checker_function=lambda x: x in TREE_ROOT_METHODS,
                equate=False,
            ),
            # root2           midlongestspan                          tree; root1 is used in
            #                minavgleafdist                          iteration 1 and 2,
            #                                                        root2 in later
            #                                                        iterations.
            _Option(
                ["-root2", "root2"],
                "Method used to root tree in iteration 2",
                checker_function=lambda x: x in TREE_ROOT_METHODS,
                equate=False,
            ),
            # seqtype         protein              auto               Sequence type.
            #                nucleo
            #                auto
            _Option(
                ["-seqtype", "seqtype"],
                "Sequence type",
                checker_function=lambda x: x in SEQUENCE_TYPES,
                equate=False,
            ),
            # smoothscoreceil Floating point       [1]                Maximum value of column
            #                                                        score for smoothing
            #                                                        purposes.
            _Option(
                ["-smoothscoreceil", "smoothscoreceil"],
                "Maximum value of column score for smoothing",
                checker_function=lambda x: isinstance(x, float),
                equate=False,
            ),
            # smoothwindow    Integer              7                  Window used for anchor
            #                                                        column smoothing.
            _Option(
                ["-smoothwindow", "smoothwindow"],
                "Window used for anchor column smoothing",
                checker_function=lambda x: isinstance(x, int),
                equate=False,
            ),
            # SUEFF           Floating point value 0.1                Constant used in UPGMB
            #                between 0 and 1.                        clustering. Determines
            #                                                        the relative fraction
            #                                                        of average linkage
            #                                                        (SUEFF) vs. nearest-
            #                                                        neighbor linkage (1
            #                                                        SUEFF).
            _Option(
                ["-sueff", "sueff"],
                "Constant used in UPGMB clustering",
                checker_function=lambda x: isinstance(x, float),
                equate=False,
            ),
            # tree1           File name            None               Save tree produced in
            _Option(["-tree1", "tree1"], "Save Newick tree from iteration 1", equate=False),
            # tree2                                                   first or second
            #                                                        iteration to given file
            #                                                        in Newick (Phylip-
            #                                                        compatible) format.
            _Option(["-tree2", "tree2"], "Save Newick tree from iteration 2", equate=False),
            # weight1         none                 clustalw           Sequence weighting
            _Option(
                ["-weight1", "weight1"],
                "Weighting scheme used in iteration 1",
                checker_function=lambda x: x in WEIGHTING_SCHEMES,
                equate=False,
            ),
            # weight2         henikoff                                scheme.
            #                henikoffpb                              weight1 is used in
            #                gsc                                     iterations 1 and 2.
            #                clustalw                                weight2 is used for
            #                threeway                                tree-dependent
            #                                                        refinement.
            #                                                        none=all sequences have
            #                                                        equal weight.
            #                                                        henikoff=Henikoff &
            #                                                        Henikoff weighting
            #                                                        scheme.
            #                                                        henikoffpb=Modified
            #                                                        Henikoff scheme as used
            #                                                        in PSI-BLAST.
            #                                                        clustalw=CLUSTALW
            #                                                        method.
            #                                                        threeway=Gotoh three-
            #                                                        way method.
            _Option(
                ["-weight2", "weight2"],
                "Weighting scheme used in iteration 2",
                checker_function=lambda x: x in WEIGHTING_SCHEMES,
                equate=False,
            ),
            # ################### FORMATS #######################################
            # Multiple formats can be specified on the command line
            # If -msf appears it will be used regardless of other formats
            # specified. If -clw appears (and not -msf), clustalw format will be
            # used regardless of other formats specified. If both -clw and
            # -clwstrict are specified -clwstrict will be used regardless of
            # other formats specified. If -fasta is specified and not -msf,
            # -clw, or clwstrict, fasta will be used. If -fasta and -html are
            # specified -fasta will be used. Only if -html is specified alone
            # will html be used. I kid ye not.
            # clw                no              Write output in CLUSTALW format (default is
            #                                   FASTA).
            _Switch(
                ["-clw", "clw"],
                "Write output in CLUSTALW format (with a MUSCLE header)",
            ),
            # clwstrict          no              Write output in CLUSTALW format with the
            #                                   "CLUSTAL W (1.81)" header rather than the
            #                                   MUSCLE version. This is useful when a post-
            #                                   processing step is picky about the file
            #                                   header.
            _Switch(
                ["-clwstrict", "clwstrict"],
                "Write output in CLUSTALW format with version 1.81 header",
            ),
            # fasta              yes             Write output in FASTA format. Alternatives
            #                                   include clw,
            #                                   clwstrict, msf and html.
            _Switch(["-fasta", "fasta"], "Write output in FASTA format"),
            # html               no              Write output in HTML format (default is
            #                                   FASTA).
            _Switch(["-html", "html"], "Write output in HTML format"),
            # msf                no              Write output in MSF format (default is
            #                                   FASTA).
            _Switch(["-msf", "msf"], "Write output in MSF format"),
            # Phylip interleaved - undocumented as of 3.7
            _Switch(["-phyi", "phyi"], "Write output in PHYLIP interleaved format"),
            # Phylip sequential - undocumented as of 3.7
            _Switch(["-phys", "phys"], "Write output in PHYLIP sequential format"),
            # ################# Additional specified output files #########
            _Option(
                ["-phyiout", "phyiout"],
                "Write PHYLIP interleaved output to specified filename",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-physout", "physout"],
                "Write PHYLIP sequential format to specified filename",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-htmlout", "htmlout"],
                "Write HTML output to specified filename",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-clwout", "clwout"],
                "Write CLUSTALW output (with MUSCLE header) to specified " "filename",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-clwstrictout", "clwstrictout"],
                "Write CLUSTALW output (with version 1.81 header) to " "specified filename",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-msfout", "msfout"],
                "Write MSF format output to specified filename",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-fastaout", "fastaout"],
                "Write FASTA format output to specified filename",
                filename=True,
                equate=False,
            ),
            # ############# END FORMATS ###################################
            # anchors            yes             Use anchor optimization in tree dependent
            #                                   refinement iterations.
            _Switch(
                ["-anchors", "anchors"],
                "Use anchor optimisation in tree dependent " "refinement iterations",
            ),
            # noanchors          no              Disable anchor optimization. Default is
            #                                   anchors.
            _Switch(
                ["-noanchors", "noanchors"],
                "Do not use anchor optimisation in tree dependent " "refinement iterations",
            ),
            # group              yes             Group similar sequences together in the
            #                                   output. This is the default. See also
            #                                   stable.
            _Switch(["-group", "group"], "Group similar sequences in output"),
            # stable             no              Preserve input order of sequences in output
            #                                   file. Default is to group sequences by
            #                                   similarity (group).
            _Switch(
                ["-stable", "stable"],
                "Do not group similar sequences in output (not supported in v3.8)",
            ),
            # ############# log-expectation profile score ######################
            # One of either -le, -sp, or -sv
            #
            # According to the doc, spn is default and the only option for
            # nucleotides: this doesnt appear to be true. -le, -sp, and -sv can
            # be used and produce numerically different logs (what is going on?)
            #
            # spn fails on proteins
            # le                 maybe           Use log-expectation profile score (VTML240).
            #                                    Alternatives are to use sp or sv. This is
            #                                    the default for amino acid sequences.
            _Switch(["-le", "le"], "Use log-expectation profile score (VTML240)"),
            # sv                 no              Use sum-of-pairs profile score (VTML240).
            #                                   Default is le.
            _Switch(["-sv", "sv"], "Use sum-of-pairs profile score (VTML240)"),
            # sp                 no              Use sum-of-pairs protein profile score
            #                                   (PAM200). Default is le.
            _Switch(["-sp", "sp"], "Use sum-of-pairs protein profile score (PAM200)"),
            # spn                maybe           Use sum-of-pairs nucleotide profile score
            #                                   (BLASTZ parameters). This is the only option
            #                                   for nucleotides, and is therefore the
            #                                   default.
            _Switch(["-spn", "spn"], "Use sum-of-pairs protein nucleotide profile score"),
            # ############# END log-expectation profile score ######################
            # quiet              no              Do not display progress messages.
            _Switch(["-quiet", "quiet"], "Use sum-of-pairs protein nucleotide profile score"),
            # refine             no              Input file is already aligned, skip first
            #                                   two iterations and begin tree dependent
            #                                   refinement.
            _Switch(["-refine", "refine"], "Only do tree dependent refinement"),
            # core               yes in muscle,  Do not catch exceptions.
            #                   no in muscled.
            _Switch(["-core", "core"], "Catch exceptions"),
            # nocore             no in muscle,   Catch exceptions and give an error message
            #                   yes in muscled. if possible.
            _Switch(["-nocore", "nocore"], "Do not catch exceptions"),
            # termgapsfull       no              Terminal gaps penalized with full penalty.
            #                                   [1] Not fully supported in this version.
            #
            # termgapshalf       yes             Terminal gaps penalized with half penalty.
            #                                   [1] Not fully supported in this version.
            #
            # termgapshalflonger no              Terminal gaps penalized with half penalty if
            #                                   gap relative to
            #                                   longer sequence, otherwise with full
            #                                   penalty.
            #                                   [1] Not fully supported in this version.
            # verbose            no              Write parameter settings and progress
            #                                   messages to log file.
            _Switch(["-verbose", "verbose"], "Write parameter settings and progress"),
            # version            no              Write version string to stdout and exit.
            _Switch(["-version", "version"], "Write version string to stdout and exit"),
        ]
        self.parameters += [
            _Option(
                ["-spscore", "spscore"],
                "Filename of multiple alignment for which to compute SP score",
                filename=True,
                equate=False,
            ),
            _Option(
                ["-scorefile", "scorefile"],
                "Write alignment scores to specified filename",
                filename=True,
                equate=False,
            ),
        ]
        super(MuscleCommandline, self).__init__(cmd, **kwargs)


class FlashCommandline(AbstractCommandline):
    """
    Interface to FLASH alignment.

    The following example:

    cline = FlashCommandline(max_overlap=200,
                read1='read1.fastq',read2='read2.fastq')
    stdout,stderr = cline()

    is equivalent to the command line:

    $ flash read1.fastq read2.fastq -M 200
    """

    def __init__(self, cmd="flash", **kwargs):
        self.program_name = cmd
        self.parameters = [
            _Argument(["read1"], "read1 input filename", filename=True),
            _Argument(["read2"], "read2 input filename", filename=True),
            _Switch(
                ["-O", "allow_outies"],
                'also try combining read pairs in the "outie" orientation',
            ),
            ## alignment
            _Option(
                ["-M", "max_overlap"],
                "max overlap length expected in 90% of read pairs",
                checker_function=lambda x: x > 0,
                equate=False,
            ),
            ## I/O, processing
            _Option(
                ["-o", "output_prefix"],
                'Prefix of output files. Default: "out".',
                equate=False,
            ),
            _Option(
                ["-d", "output_directory"],
                "Path to directory for output files. Default: current working directory.",
                equate=False,
            ),
            _Option(
                ["-t", "threads"],
                "Set the number of worker threads. This is in addition to the I/O threads.",
                equate=False,
            ),
        ]
        super(FlashCommandline, self).__init__(cmd, **kwargs)
