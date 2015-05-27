***  MSplazer - Generic Multi-Split Chaining ***
http://www.seqan.de/projects/MSplazer.html

---------------------------------------------------------------------------
Table of Contents
---------------------------------------------------------------------------
  1.   Overview
  2.   Installation
  3.   Usage
  4.   Output Format
  5.   Example
  6.   Contact

---------------------------------------------------------------------------
1. Overview
---------------------------------------------------------------------------


MSplazer (MultiSplitRazerS) is a tool primarily designed for split mapping of 
sequencing reads.

MSplazer is a sound generic multi-split chaining method using the
C++ library SeqAn that uses SeqAns exact local aligner Stellar to detect splits of a
read. Compatible local matches of a read are then identified, and all compatibility
information is stored in a split-read graph representation of the matches. We then use a
DAG shortest path algorithm to determine the most probable chain of splits, and report
the underlying breakpoints.



---------------------------------------------------------------------------
2. Installation
---------------------------------------------------------------------------

A precompiled linux 64-bit binary of MSplazer can be downloaded from the 
SeqAn projects download page: http://www.seqan.de/downloads/projects.html

MSplazer is distributed with SeqAn - The C++ Sequence Analysis Library (see 
http://www.seqan.de). To build SplazerS yourself you can check out the latest 
SVN version of MSplazer and SeqAn with:

  1)  svn co http://svn.mi.fu-berlin.de/seqan/trunk/seqan seqan-trunk
  2)  mkdir seqan-trunk/build/
  3)  cd seqan-trunk/build/
  4)  cmake .. -DCMAKE_BUILD_TYPE=Release
  5)  make msplazer
  6)  ./sandbox/ktrappe/apps/msplazer/msplazer --help

Alternatively, you can download the latest SeqAn snapshot and build MSplazer
with:

  1)  Download the latest snapshot of SeqAn
  2)  Unzip it to a directory of your choice (e.g. snapshot)
  3)  cd snapshot/apps
  4)  make msplazer

If the make command was successful, you can skip ahead to step 5. Otherwise,
you may have to manually add MSplazer to the Makefile targets. For this, open 
the file "Makefile" and add "msplazer/msplazer" to the "TARGETS = " line. 
Also, add a line "msplazer: check_seqan_base msplazer/msplazer" anywhere below
the "TARGETS = " and above the "check_seqan_base:" line. Now the make command 
should be successful.

  5)  cd msplazer
  6)  ./msplazer --help

If succesful, an executable file msplazer was built and a brief usage 
description is dumped.


---------------------------------------------------------------------------
3. Usage
---------------------------------------------------------------------------

To get a short usage description of SplazerS, you can execute msplazer -h or 
msplazer --help.

Usage: msplazer -d <FASTA sequence file> -q <FASTA sequence file> [Options]
       msplazer -d <FASTA sequence file> -q <FASTA sequence file> -m <GFF match file> [Options]

  Using the first line, MSplazer will first run Stellar internally on the given
  input files. Using the --matchfile (-m) option, MSplazer skips this step and
  directly uses the given matches, preferable precalculated using Stellar, from
  the GFF file.

  To get a short usage description of MSplazer, you can execute msplazer -h or 
  msplazer --help.

---------------------------------------------------------------------------
3.1. Non-optional arguments
---------------------------------------------------------------------------

  MSplazer, like STELLAR, always expects the two parameters --database and --query (or -d and
  -q) to be set.

  [ -d FILE ],  [ --database FILE ]
 
  Set the name of the file containing the database sequence(s) in
  (multi-)Fasta format.  
 
  [ -q FILE ],  [ --query FILE ]
 
  Set the name of the file containing the query sequence(s) in (multi-)
  Fasta format. All query sequences will be compared to all database
  sequences. Important: Following conventions, all Ids (line starting 
  with '>') have to be unique already within the first part until the 
  first whitespace! For example:
  >MSplazer|group=6|reads=9 readId=1
  >MSplazer|group=6|reads=9 readId=2
  are not unique. 

  Without any additional parameters, MSplazer would call Stellar and then chain
  those matches of each read, that have either an overlap of at least 0.5 (50% 
  of each match length) or a distance between the matches of at most 10 bp, 
  and that miss at most 15 bp at the end or beginning of the read.
  The program calls Stellar with default options: Stellar would 
  compare the query sequence(s) to both strands of the database sequence(s) 
  with an error rate of 0.05 (i.e. 5% errors, an identity of 95%) and a 
  minimal length of 100, and dump all local alignments in an output file named
  "stellar.gff" (see Stellar Readme).

  Two output files will be generated: breakpoint.gff includes all structural 
  variant breakpoints in GFF format, indel.gff includes small indel (usually 
  <10bp) in GFF format.

  The default behaviour can be modified by adding the following options to
  the command line.

---------------------------------------------------------------------------
3.2. Main Options
---------------------------------------------------------------------------

---------------------------------------------------------------------------
3.2.1. Stellar-Inherited Main Options
---------------------------------------------------------------------------

 [ -e REAL ],  [ --epsilon REAL ]
  
  Set the maximal error rate for local alignments. REAL must be a floating
  point number between 0 and 0.25. The default value is 0.05. REAL is the
  number of edit operations needed to transform the aligned query substring
  to the aligned database substring divided by the length of the local
  alignment. For example specify '-e 0.1' for a maximal error rate of 10%
  (90% identity).

  [ -l NUM ],  [ --minLength NUM ]
  
  Set the minimal length of a local alignment. The default value is 100.
  NUM is the minimal length of an alignment, not the length of the
  aligned substrings.

  [ -f ],  [ --forward ]

  Only compare query sequence(s) to the positive/forward strand of the
  database sequence(s). By default, both strands are scanned.

  [ -r ],  [ --reverse ]

  Only compare query sequence(s) to the negative/reverse-complemented 
  database sequence(s). By default, both strands are scanned.
 
  [ -v ],  [ --verbose ]
  
  Verbose. Print extra information and running times.
  For more Stellar options see Stellar documentation/readme.

---------------------------------------------------------------------------
3.2.2. MSplazer-Specific Main Options
---------------------------------------------------------------------------

  [ -m FILE ],  [ --matchfile FILE ]
 
  Set the name of a file containing matches in GFF format. 

  When setting option --matchfile (-m), MSplazer uses the matches from the 
  given input file instead of calling Stellar (practical for running multiple
  tests on the same data). See also sample file stellar.gff.
 
  [ -oth REAL ],  [ --overlapThresh REAL ]
  
  Required overlap for matches of one read to be considered for chaining 
  (main criterion, default 0.5).

  [ -gth NUM ],  [ --gapThresh NUM ]
  
  Maximal allowed distance between two matches to be considered for 
  chaining. This criterion only applies if the matches do not overlap 
  (default 10).

  [ -ith NUM], [ --initGapThresh NUM ]
 
  Maximal allowed length of leading or ending gap for the whole read(!) 
  (default 15).

  [ -st NUM], [ --support NUM ]
  
  Required number of supporting reads (or contigs) for a breakpoint to be
  written to the output file (see breakpoint output, default 2).

  [ -tp NUM], [ --transPen NUM ]
 
  Interchromosomal translocation penalty (default 5). Added to edge weight
  between matches that map to different database sequences.

  [ -ip NUM], [ --invPen NUM ]
 
  Inversion penalty (default 5). Added to edge weight between matches that 
  map to different strands of the same database sequence. Note: Penalty is
  only added if there is no transPen already on the edge.

  [ -op NUM], [ --orderPen NUM ]

  Order penalty (default 5). Added to edge weight between matches that 
  map to a database sequence in a different order than they are in the read 
  sequence. Note: Penalty is  only added if there is no transPen  or invPen 
  already on the edge.



---------------------------------------------------------------------------
3.3 Output Options
---------------------------------------------------------------------------

  [ -i ],  [ --outdir ]

  Optional existing(!) output directory for all output files (breakpoint.gff,
  indel.gff, *.dot).

  [ -j ],  [ --jobName ]
  
  Optional jobname that will be added to all output files in the format
  jobname_breakpoint.gff
  jobname_indel.gff
  read#_jobname.dot. # stands for the sequence number of the read within the 
  query input file.
  
  [ -do ], [ --dots ]
  
  Switches output of DOT files on and off (default off). Each DOT file
  includes the graph representation of one read/query.

---------------------------------------------------------------------------
4. Output Formats
---------------------------------------------------------------------------

MSplazer currently supports the GFF output format for reporting brakpoints. 
 
---------------------------------------------------------------------------
4.1. General Feature Format (GFF)
---------------------------------------------------------------------------

The General Feature Format is specified by the Sanger Institute as a tab-
delimited text format with the following columns:

<seqname> <src> <feat> <start> <end> <score> <strand> <frame> [attr] [cmts]

See also: http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml
Consistent with this specification MSplazer GFF output looks as follows:

GNAME msplazer SVTYPE GBEGIN GEND . GSTRAND . ATTRIBUTES
StartSeqId      Label   SV type sPos    sPos    .       +/-     .       
Match value description:

  GNAME        Name of the genome sequence (see --genome-naming)
  msplazer     Constant
  SVTYPE       SV type (deletion, insertion, inversion, translocation)
  GBEGIN       Beginning position in the genome sequence 
               (positions are counted from 1)
  GEND         End position in the genome sequence (included!)
  .            Constant
  GSTRAND      '+'=forward strand or '-'=reverse strand
  .            Constant
  ATTRIBUTES   A list of attributes in the format <tag_name>[=<tag>] 
               separated by ';'

Attributes are: 

  ID=          Random number as identifier
  size=        Indel size (deletion and insertion only), "-" before the 
               size indicates an insertion
  seq=         Sequence content of insertion (obviously insertion only).
  endChr=      Name of the second genome sequence (inversion and 
               interchromosomal translocation only)
  endPos=      First position in second genome sequence (inversion and 
               interchromosomal translocation only)
  endStrand=   '+'=forward strand or '-'=reverse strand of second genome 
               sequence (inversion and interchromosomal translocation only)
  support=     Number of supporting reads for this breakpoint 
  supportIds=  IDs of supporting reads

For matches on the reverse strand, GBEGIN and GEND are positions on the 
related forward strand. It holds GBEGIN < GEND, regardless of GSTRAND.


---------------------------------------------------------------------------
5. Example
---------------------------------------------------------------------------

The subfolder 'example' contains a little example. The genome file is chr10.fa,
the sample read (or better contig) file 20_contigs_chr10_sim_uniqIds.fa. There
is also a file 20_contigs_chr10_sim_uniqIds.gff that contains pre-calculated
Stellar matches for the mentioned input files.

The default calls would then be

../msplazer -d chr10.fa -q 20_contigs_chr10_sim_uniqIds.fa
or
../msplazer -d chr10.fa -q 20_contigs_chr10_sim_uniqIds.fa -m 20_contigs_chr10_sim_uniqIds.gff

Both calls produce the default output files breakpoint.gff and indel.gff. 
Note that both files should be empty since the default support is 2. The
file generate_outputs.sh contains the calls for Stellar and MSplazer with
some specified optional parameters that directs the output to the directories
internalStellar and externalStellar.

---------------------------------------------------------------------------
6. Contact
---------------------------------------------------------------------------

For questions or comments, contact:
  Kathrin Trappe  <kathrin.trappe@fu-berlin.de>
