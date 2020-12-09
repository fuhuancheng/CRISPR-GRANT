Introduction
============

GRACIOUS, a stand-alone graphical CRISPR indel analysis tool, could be
easily installed for multi-platform including Linux, Windows, and MacOS.
GRACIOUS offered a straightforward GUI by simple click-and-run for
genome editing analysis of single or pooled amplicons and one-step
analysis for whole-genome sequencing. Moreover, it also exhibited
shorter run-time compared with tools currently available. Therefore,
GRACIOUS is a valuable addition to the current CRISPR toolkits that
significantly lower the barrier for wet-lab researchers to conduct indel
analysis from large NGS datasets.

GRACIOUS requires input from raw FASTQ sequencing data and reference
sequence in FASTA format. Fastp will be first used for quality check of
the input data. After that, qualified reads will be mapped to the
reference sequence by BWA-MEM. The resulting SAM file was then converted
to BAM by samtools, which will be subsequently used for consensus and
variants analysis by VarScan2. The output data contains QC reports,
FASTQ of qualified reads, reference-mapping results, table of consensus
and variants, and publication-quality plots for the number of modified
and unmodified reads, alignment of top numbered reads to reference, and
the frequency of indels at each position.

Download and install
====================

System requirements:

-   Operating systems: Windows 7 or 10; macOS 10.13 or later; GNU/Linux
    (Ubuntu 16.04 or later). The program had been tested on Windows 7,
    Windows 10, macOS 10.13 (High Sierra), Ubuntu 16.04 and deepin
    linux 20.

-   CPU: x86_64 CPUs (64 bit).

-   RAM: at least 4 GB RAM, 8 GB is recommended. For very large or whole
    genome sequencing (WGS) dataset, larger RAM may be necessary.

-   Disk space: it depends on the sequencing data size. For WGS data of
    mouse (*Mus musculus*), as large as 500GB disk space would be
    required.

Download the GRACIOUS package from website
(<https://github.com/fuhuancheng/GRACIOUS/releases>), uncompressed it,
go to the folder and click to open the program. The GUI would look like
Fig[1](#fig:mac-gui).

![GUI example](doc/figures/mac-gui2.png)

The usage of GRACIOUS is simple and intuitive.

Usage for amplicon or pooled amplicon analysis
==============================================

Required files or inputs for analysis:

-   FASTQ file(s) from single-end or paired-end sequencing.

-   Reference sequence in FASTA file format.

-   Output folder for analysis results.

Analysis for amplicon or pooled amplicons includes steps
(Fig[2](#fig:amplicon-gui)):

1.  Press \"Browse\" button to select FASTQ file 1 and 2.

2.  Press \"Browse\" button to select reference file in FASTA format.

3.  Press \"Browse\" button to select folder for result files.

4.  Choose additional parameters if needed, or just leave it as default.

5.  Press \"begin analysis\" button to begin analysis.

6.  When finished, go to the result folder to check the result files.

![Steps for amplicon(s)
analysis](doc/figures/amplicon-gui.png)

Usage for WGS (whole genome sequencing) analysis
================================================

Required inputs for analysis:

-   FASTQ file(s) from single-end or paired-end sequencing.

-   Reference sequence in FASTA file format.

-   Output folder for analysis results.

-   A file containing specific regions.

Since the reference genomes are generally huge, the genome reference
files were not provided with the installer packages. However, it's easy
to download genome sequence for a variety of organisms from UCSC
(<https://hgdownload.soe.ucsc.edu/downloads.html>) or NCBI
(<https://www.ncbi.nlm.nih.gov/genome/>). The downloaded genome sequence
should be in FASTA format (see
section[4.1.2](#sec:fasta)).

Analysis for WGS include steps
(Fig[3](#fig:wgs-gui)):

1.  Press \"Browse\" button to select FASTQ file 1 and 2.

2.  Press \"Browse\" button to select reference file in FASTA format.

3.  Press \"Browse\" button to select folder for result files.

4.  Choose additional parameters if needed, or just leave it as default.

5.  Press \"begin analysis\" button to begin analysis.

6.  When finished, go to the result folder to check the result files.

![Steps for WGS analysis](doc/figures/wgs-gui.png)

Input file format
-----------------

### FASTQ: sequencing file format

FASTQ files (\*.fastq) or compressed FASTQ files (\*.fastq.gz) were both
accepted.

### FASTA: reference file format {#sec:fasta}

The reference file must be in FASTA file format, file extension is
arbitrary. The detailed FASTA file format explaination could be found
from <https://en.wikipedia.org/wiki/Fasta_format>.

A sequence in FASTA format begins with a single-line description,
followed by lines of sequence data. The definition line (defline) is
distinguished from the sequence data by a greater-than (\>) symbol at
the beginning. The word following the \"\>\" symbol is the identifier of
the sequence, and the rest of the line is the description (optional).
There should be **NO** space between the \"\>\" and the first letter of
the identifier. It is recommended that all lines of text be shorter than
or equal to 80 characters in length. Below is an example sequence in
FASTA format:

    >XM_030254057.1 PREDICTED: Mus musculus actin, beta (Actb)
    ATGCACTGAGGCACGGCAGGCCCAGAGCATCTCACCTGAAGCACCCTTCTTGCCTAAATCCAGCTTTCTG
    TCACACTCTCCCAGAAGGAGGGGAGAGGGGGTAAAAAAATGCTGCACTGTGCGGCGAGGCCGGTGAGTGA
    GCGACGCGGAGCCAATCAGCGCCCGCCGTTCCGAAAGTTGCCTTTTATGGCTCGAGTGGCCGCTGTGGCG
    TCCTATAAAACCCGGCGGCGCAACGC...

In cases such as reference for amplicon pools, all the sequence names in
the given reference file should **NOT** be in duplicate.

Result files
------------

### Alignment of top counts

In the result folder, the alignment of top count sequences will be
given, as shown in fig[4](#fig:align).

![Example of alignment plot](doc/figures/top_alignment.png)

### Indel frequency along sequence

In the result folder, the plot of indel frequency along sequence will
also be given, as shown in
figure[5](#fig:indel-freq).

![Example of indel frequency along
sequence](doc/figures/freq-plot.png)

### Basic statistics

Performance
===========

The performance of GRACIOUS and other relating software was tested on
Ubuntu 16.04 (CPU: Intel i7-8700K).

Command line usage
==================

GRACIOUS also provides command line usage for some situations when
needed.

    indel_analysis -1=FASTQ1/file/path -2=FASTQ2/file/path -r=reference/file/path -o=output/file/path -t=CPU_cores

Command \"indel_analysis -h\" would give more detailed help usage
information.

Compilation
===========

Requirements:

-   gcc or clang compiler. On Windows, Mingw-w64 should be installed and
    used for compilation.

-   Nim (greater than Version 1.2.0). Please consult
    <https://nim-lang.org> for installation.

-   ggplotnim (greater than version 0.3.18). ggplotnim could be
    installed via `nimble install ggplotnim`.

-   ui (greater than version 0.9.4). ui could be installed via
    `nimble install ui`.

Download the source code, uncompressed it and change into the source
code directory.

On Linux:
      make unix
      make linux

On Mac:
      make unix
      make mac

On Windows:
      make windows
      make win

Then compressed GRACIOUS package in zip format will then be in bin/ folder.

You also need fastp, flash, bwa, samtools, VarScan2 installed in bin
folder.

The resulting bin folder should be like below:

      bin/
      |-- bwa/
      |-- flash/
      |-- jre11/
      |-- mafft/
      |-- samtools/
      |-- varscan2/
      |-- csv_to_fasta
      |-- fasta_to_plot
      |-- fastp
      |-- indel_analysis
      |-- reads_count_plot
      |-- sam_count
      |-- varToCsv
      |-- var_plot
      |-- wgsSubRegion

Get help
========

If you have any questions or encountered bugs when using GRACIOUS,
please feel free to contact us by email
([fuhuancheng\@foxmail.com](fuhuancheng@foxmail.com)).

Additional parameters for FLASH
===============================

Additional options, if necessary, can be added and passed directly to
internally used tools in Additional Options tab
(Fig[6](#fig:additional-options)).

![GUI of additional
options](doc/figures/additional_options.png)

Here are all available additional options for flash. Detailed additional
parameters for FLASH can be accessed by execute command \"flash -h\" in
terminal.

    -m, --min-overlap=NUM   The minimum required overlap length between two
                         reads to provide a confident overlap.  Default:
                         10bp.

    -M, --max-overlap=NUM   Maximum overlap length expected in approximately
                         90% of read pairs.  It is by default set to 65bp,
                         which works well for 100bp reads generated from a
                         180bp library, assuming a normal distribution of
                         fragment lengths.  Overlaps longer than the maximum
                         overlap parameter are still considered as good
                         overlaps, but the mismatch density (explained below)
                         is calculated over the first max_overlap bases in
                         the overlapped region rather than the entire
                         overlap.  Default: 65bp, or calculated from the
                         specified read length, fragment length, and fragment
                         length standard deviation.

    -x, --max-mismatch-density=NUM
                         Maximum allowed ratio between the number of
                         mismatched base pairs and the overlap length.
                         Two reads will not be combined with a given overlap
                         if that overlap results in a mismatched base density
                         higher than this value.  Note: Any occurence of an
                         'N' in either read is ignored and not counted
                         towards the mismatches or overlap length.  Our
                         experimental results suggest that higher values of
                         the maximum mismatch density yield larger
                         numbers of correctly merged read pairs but at
                         the expense of higher numbers of incorrectly
                         merged read pairs.  Default: 0.25.

    -O, --allow-outies      Also try combining read pairs in the "outie"
                         orientation, e.g.

                              Read 1: <-----------
                              Read 2:       ------------>

                         as opposed to only the "innie" orientation, e.g.

                              Read 1:       <------------
                              Read 2: ----------->

                         FLASH uses the same parameters when trying each
                         orientation.  If a read pair can be combined in
                         both "innie" and "outie" orientations, the
                         better-fitting one will be chosen using the same
                         scoring algorithm that FLASH normally uses.

                         This option also causes extra .innie and .outie
                         histogram files to be produced.

    -p, --phred-offset=OFFSET
                         The smallest ASCII value of the characters used to
                         represent quality values of bases in FASTQ files.
                         It should be set to either 33, which corresponds
                         to the later Illumina platforms and Sanger
                         platforms, or 64, which corresponds to the
                         earlier Illumina platforms.  Default: 33.

    -r, --read-len=LEN
    -f, --fragment-len=LEN
    -s, --fragment-len-stddev=LEN
                         Average read length, fragment length, and fragment
                         standard deviation.  These are convenience parameters
                         only, as they are only used for calculating the
                         maximum overlap (--max-overlap) parameter.
                         The maximum overlap is calculated as the overlap of
                         average-length reads from an average-size fragment
                         plus 2.5 times the fragment length standard
                         deviation.  The default values are -r 100, -f 180,
                         and -s 18, so this works out to a maximum overlap of
                         65 bp.  If --max-overlap is specified, then the
                         specified value overrides the calculated value.

                         If you do not know the standard deviation of the
                         fragment library, you can probably assume that the
                         standard deviation is 10% of the average fragment
                         length.

    --cap-mismatch-quals    Cap quality scores assigned at mismatch locations
                         to 2.  This was the default behavior in FLASH v1.2.7
                         and earlier.  Later versions will instead calculate
                         such scores as max(|q1 - q2|, 2); that is, the
                         absolute value of the difference in quality scores,
                         but at least 2.  Essentially, the new behavior
                         prevents a low quality base call that is likely a
                         sequencing error from significantly bringing down
                         the quality of a high quality, likely correct base
                         call.

    -z, --compress          Compress the output files directly with zlib,
                         using the gzip container format.  Similar to
                         specifying --compress-prog=gzip and --suffix=gz,
                         but may be slightly faster.

    --compress-prog=PROG    Pipe the output through the compression program
                         PROG, which will be called as `PROG -c -',
                         plus any arguments specified by --compress-prog-args.
                         PROG must read uncompressed data from standard input
                         and write compressed data to standard output when
                         invoked as noted above.
                         Examples: gzip, bzip2, xz, pigz.

    --compress-prog-args=ARGS
                         A string of additional arguments that will be passed
                         to the compression program if one is specified with
                         --compress-prog=PROG.  (The arguments '-c -' are
                         still passed in addition to explicitly specified
                         arguments.)

Additional parameters for BWA
=============================

Additional parameters for FLASH can be accessed by execute command \"bwa
mem\" in terminal. Here are all available additional options for bwa
mapping:

            -k INT        minimum seed length [19]
            -w INT        band width for banded alignment [100]
            -d INT        off-diagonal X-dropoff [100]
            -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]
            -y INT        seed occurrence for the 3rd round seeding [20]
            -c INT        skip seeds with more than INT occurrences [500]
            -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [0.50]
            -W INT        discard a chain if seeded bases shorter than INT [0]
            -m INT        perform at most INT rounds of mate rescues for each read [50]
            -S            skip mate rescue
            -P            skip pairing; mate rescue performed unless -S also in use

    Scoring options:

            -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [1]
            -B INT        penalty for a mismatch [4]
            -O INT[,INT]  gap open penalties for deletions and insertions [6,6]
            -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
            -L INT[,INT]  penalty for 5'- and 3'-end clipping [5,5]
            -U INT        penalty for an unpaired read pair [17]

            -x STR        read type. Setting -x changes multiple parameters unless overridden [null]
                          pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)
                          ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)
                          intractg: -B9 -O16 -L5  (intra-species contigs to ref)

    Input/output options:

            -p            smart pairing (ignoring in2.fq)
            -R STR        read group header line such as '@RG\tID:foo\tSM:bar' [null]
            -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]
            -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)
            -5            for split alignment, take the alignment with the smallest coordinate as primary
            -q            don't modify mapQ of supplementary alignments
            -K INT        process INT input bases in each batch regardless of nThreads (for reproducibility) []
            -v INT        verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [3]
            -T INT        minimum score to output [30]
            -h INT[,INT]  if there are <INT hits with score >80% of the max score, output all in XA [5,200]
            -a            output all alignments for SE or unpaired PE
            -C            append FASTA/FASTQ comment to SAM output
            -V            output the reference FASTA header in the XR tag
            -Y            use soft clipping for supplementary alignments
            -M            mark shorter split hits as secondary
            -I FLOAT[,FLOAT[,INT[,INT]]]
                          specify the mean, standard deviation (10% of the mean if absent), max
                          (4 sigma from the mean if absent) and min of the insert size distribution.
                          FR orientation only. [inferred]

Additional parameters for VarScan2
==================================

Additional parameters for FLASH can be accessed by execute command
\"java -jar VarScan.jar\" in terminal. Here are all available additional
options for calling consensus and variants:

    --min-coverage  Minimum read depth at a position to make a call [8]
    --min-reads2    Minimum supporting reads at a position to call variants [2]
    --min-avg-qual  Minimum base quality at a position to count a read [15]
    --min-var-freq  Minimum variant allele frequency threshold [0.01]
    --min-freq-for-hom  Minimum frequency to call homozygote [0.75]
    --p-value   Default p-value threshold for calling variants [99e-02]
    --strand-filter Ignore variants with >90% support on one strand [1]
    --variants  Report only variant (SNP/indel) positions [0]
