#include <seqan/seq_io.h>

seqan::FaiIndex faiIndex;
int res = read(faiIndex, "path/to/file.fasta");
if (res != 0)
    std::cerr << "ERROR: Could not read FAI index path/to/file.fasta.fai\n";
#include <seqan/seq_io.h>

seqan::FaiIndex faiIndex;
int res = read(faiIndex, "path/to/file.fasta", "path/to/index.fai");
if (res != 0)
    std::cerr << "ERROR: Could not load FAI index path/to/index.fai\n";

unsigned idx = 0;
if (!getIdByName(faiIndex, "chr1", idx))
    std::cerr << "ERROR: FAI index has no entry for chr1.\n";

unsigned seqLength = sequenceLength(faiIndex, idx);

// Load first 1000 characters of chr1.
seqan::CharString seqChr1Prefix;
if (readRegion(seqChr1Prefix, faiIdx, idx, 0, 1000) != 0)
    std::cerr << "ERROR: Could not load chr1.\n";

// Load all of chr1.
seqan::CharString seqChr1;
if (readSequence(seqChr1, faiIdx, idx) != 0)
    std::cerr << "ERROR: Could not load chr1.\n";

