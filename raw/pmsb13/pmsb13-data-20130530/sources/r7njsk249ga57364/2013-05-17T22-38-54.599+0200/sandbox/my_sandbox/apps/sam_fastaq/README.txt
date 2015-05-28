Flexbar — flexible barcode and adapter removal, version 2.33
Bioinformatics in Quantitative Biology at BIMSB, GPLv3

=== Installation ===

To run binaries, make sure that the TBB library (Intel Threading Building Blocks) is available to the system. Flexbar binaries are provided for Linux 64, Mac OSX ≥ 10.7.4, Windows 32 and 64 bit systems on: flexbar.sourceforge.net

### Linux
One possibility is to put the file libtbb.so.2 in your working directory. To use it permanently, install the library runtime files from the repository of your distribution, or copy libtbb.so.2 to the shared library directory with the following command, or ask the administrator to install it:
sudo cp FLEXBAR_DIR/libtbb.so.2 /usr/local/lib

Or adjust the lib search path to include the directory of the lib file for the current terminal session:
export LD_LIBRARY_PATH=FLEXBAR_DIR

### Mac OSX
It applies the same as for Linux. Make the file libtbb.dylib available:
sudo cp FLEXBAR_DIR/libtbb.dylib /usr/local/lib

Or set the lib search path accordingly:
export DYLD_LIBRARY_PATH=FLEXBAR_DIR

### Windows
Keep the file tbb.dll in the directory of the Flexbar executable. Visual Studio 10 sp1 has to be installed. Alternatively, those who have not this program version can download the Visual Studio 10 sp1 redistributable package from Microsoft.
Win 32: www.microsoft.com/en-us/download/details.aspx?id=8328
Win 64: www.microsoft.com/en-us/download/details.aspx?id=13523


=== Program usage ===

Flexbar needs at least one file with sequencing reads in fasta/q or csfasta/q format as input. Further, the target name and format of reads have to be specified. For barcode based read seperation and adapter removal, a file in fasta format with barcode or adapter sequences should be provided additionally.

Please refer to the help screen (flexbar -h) or documentation on:
sourceforge.net/p/flexbar/wiki

SYNOPSIS
    flexbar -t target -f format -r reads [-b barcodes] [-a adapters] [options]

DESCRIPTION
    -h, --help
          Displays this help message.
        --version
          Display version information.
    -c, --cite
          Shows citation information.

  Basic options:
    -n, --threads NUM
          Number of threads. Default: 1.
    -t, --target STR
          Prefix for output file names.
    -r, --reads FILE
          Input file with reads, that may contain barcodes.
    -p, --reads2 FILE
          Second input file for paired read scenario.
    -f, --format STR
          Input format of reads: csfasta, csfastq, fasta, fastq, fastq-sanger,
          fastq-solexa, fastq-i1.3, fastq-i1.5, fastq-i1.8 (illumina 1.8+).

  Barcode detection:
    -b, --barcodes FILE
           Fasta file with barcodes, specify (br) to use separate barcode reads.
    -br, --barcode-reads FILE
           Fasta or fastq file with separate barcode reads, if not within reads.
    -be, --barcode-trim-end STR
           Type of barcoding within reads, see section trim-end modes. Default: ANY.
    -bn, --barcode-tail-length NUM
           Number of bases for tail trim-end modes. Default: barcode length.
    -bo, --barcode-min-overlap NUM
           Minimum overlap of barcode and read. Default: barcode length.
    -bt, --barcode-threshold NUM
           Allowed mismatches and indels per 10 bases for barcode. Default: 1.0.
    -bk, --barcode-keep
           Keep barcodes within reads instead of removal.
    -bm, --barcode-match NUM
           Alignment match score. Default: 1.
    -bi, --barcode-mismatch NUM
           Alignment mismatch score. Default: -1.
    -bg, --barcode-gap NUM
           Alignment gap score. Default: -7.

  Adapter removal:
    -a, --adapters FILE
           Fasta file with adapters, or barcodes to remove within reads.
    -as, --adapter-seq STR
           Single adapter sequence as alternative to adapters option.
    -ae, --adapter-trim-end STR
           Type of adapter removal, see section trim-end modes. Default: RIGHT.
    -an, --adapter-tail-length NUM
           Number of bases for tail trim-end modes. Default: adapter length.
    -ao, --adapter-min-overlap NUM
           Minimum overlap of adapter and read. Default: 1.
    -at, --adapter-threshold NUM
           Allowed mismatches and indels per 10 bases for adapter. Default: 3.0.
    -am, --adapter-match NUM
           Alignment match score. Default: 1.
    -ai, --adapter-mismatch NUM
           Alignment mismatch score. Default: -1.
    -ag, --adapter-gap NUM
           Alignment gap score. Default: -7.

  Filtering and trimming:
    -u, --max-uncalled NUM
          Allowed uncalled bases (N or .) in reads. Default: 0.
    -x, --pre-trim-left NUM
          Trim specified number of bases on 5' end of reads before detection.
    -y, --pre-trim-right NUM
          Trim specified number of bases on 3' end of reads before detection.
    -q, --pre-trim-phred NUM
          Trim 3' end until specified or higher quality reached.
    -z, --post-trim-length NUM
          Trim to specified read length from 3' end after removal.
    -m, --min-readlength NUM
          Minimum read length to remain after removal. Default: 18.

  Logging and tagging:
    -l, --log-level STR
          Print valid sequence alignments of reads. One of ALL, MOD, and TAB.
    -s, --single-reads
          Output single reads for partially too short paired reads.
    -o, --fasta-output
          Prefer non-quality formats fasta and csfasta for output.
    -i, --length-dist
          Write length distribution for read output files.
    -e, --number-tags
          Replace read tags by ascending number to save space.
    -g, --removal-tags
          Tag reads that are subject to adapter or barcode removal.
    -d, --random-tags
          Random read tags at barcode or adapter positions with N.

  Trim-end modes:
    ANY: longer part of read remains
    LEFT: align <= read end, right part remains
    RIGHT: align >= read start, left part remains
    LEFT_TAIL: consider first n bases of reads in alignment
    RIGHT_TAIL: use only last n bases, see tail-length options

EXAMPLES
    flexbar -t target -f fastq-i1.3 -r reads.fastq -b bar.fasta -a adap.fasta
    flexbar -t target -f csfastq -r reads.csfastq -a adapters.fasta -ae LEFT

In the first example, barcoded reads in illumina version 1.3 fastq format are demultiplexed by specifying a file with barcodes in fasta format. After read seperation based on barcodes, adapters given in fasta format are removed. In the second example, adapters in fasta format are removed if aligning before end of color-space reads, that have quality scores (csfastq). After removal the right side of reads is kept. Remaining reads are written to the file target.csfastq in same format.

Although default parameters of Flexbar are optimized to deliver good results in a large number of scenarios, the adjustment of parameters might improve results, e.g. --adapter-min-overlap and --adapter-threshold.


=== Building from SVN ===

It should be possible to compile Flexbar for any platform:

1) Check out the SVN repository to a local directory FLEXBAR_DIR, e.g. via the following command:
svn co svn://svn.code.sf.net/p/theflexibleadap/code/trunk FLEXBAR_DIR

2) Download TBB library, if you dont have Linux, Windows or Mac OSX running. For these systems the lib files are supplied together with binaries. Download the latest stable source release. It should work with version >= 3.0, then unpack the archive and run gmake in the unpacked folder.
http://www.threadingbuildingblocks.org/file.php?fid=77

3) Make the TBB library available in your library searchpath:
- For Linux 64, Mac OSX or Windows follow the steps for binaries above.
- If you compiled TBB yourself copy the compiled lib to your library searchpath and change the line "LINK_DIRECTORIES(${FLEXBAR_SOURCE_DIR}/lib/linux64)" in file FLEXBAR_DIR/src/CMakeLists.txt to include your TBB_INSTALL_DIR/build/release folder.

4) Get cmake from cmake.org and install it. Change to FLEXBAR_DIR on command line and type the following including the dot:
cmake .

CMake also allows you to create make and project files for development evironments via the -G switch, e.g.:
cmake -G "Eclipse CDT4 - Unix Makefiles" .

5) Compile the source code by issuing make in FLEXBAR_DIR. In general, the seqan and tbb library (in FLEXBAR_DIR/lib) need to be available to the compiler and linker. In case of eclipse, import the project from FLEXBAR_DIR and compile Flexbar after setting the lib path in the project settings.


=== Project folders ===

lib:      shared tbb libs for different platforms (Linux 64, Mac OSX, Windows)
include:  adapted versions of SeqAn and tbb libraries
test:     small test datasets for testing Flexbar after modifications

To run Flexbar with the test dataset for verification, make sure flexbar is reachable via the path variable and run flexbar_validate.sh within the test folder.
