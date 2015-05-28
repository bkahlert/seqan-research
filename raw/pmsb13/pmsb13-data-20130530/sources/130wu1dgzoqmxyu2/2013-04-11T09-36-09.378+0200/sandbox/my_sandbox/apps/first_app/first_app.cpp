#include <seqan/seq_io.h>

seqan::FaiIndex faiIndex;
int res = read(faiIndex, "path/to/file.fasta");
if (res != 0)
    std::cerr << "ERROR: Could not read FAI index path/to/file.fasta.fai\n";
