#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>


int main(int argc, char *argv[])
{
    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " <fasta-file> <id of sequence> <begin position> <end position>" << std::endl;
        return 1;
    }

    seqan::FaiIndex faiIndex;
    int res = read(faiIndex, argv[1]);
    if (res != 0)
    {
        std::cerr << "ERROR: Could not read FAI index path/to/file.fasta.fai\n";
        return 1;
    }
    unsigned idx = 0;
    if (!getIdByName(faiIndex, argv[2], idx))
    {
        std::cerr << "ERROR: FAI index has no entry for chr1.\n";
        return 1;
    }
    
    unsigned seqLength = sequenceLength(faiIndex, idx);

    // Load chosen characters of chromosome.
    seqan::CharString seqChr;
    if (readRegion(seqChr, faiIdx, idx, argv[3], argv[4]) != 0)
        std::cerr << "ERROR: Could not load chr1.\n";

    std::cout << "Sequence: " << seqChr << std::endl;

    exit 0;
    
}
