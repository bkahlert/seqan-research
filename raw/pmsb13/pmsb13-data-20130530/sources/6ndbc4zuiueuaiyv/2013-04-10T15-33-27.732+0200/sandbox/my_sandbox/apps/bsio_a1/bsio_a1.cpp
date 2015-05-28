#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>


int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <fasta-file>" << std::endl;
        return 1;
    }
    seqan::CharString id;
    seqan::Dna5String seq;
    seqan::CharString quals;
    seqan::SequenceStream seqStream(argv[1]);
    
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }
    
    while (!atEnd(seqStream))
    {
        // Read next sequence from the file.  In case of sequences with qualities,
        // the qualities are directly stored in the Dna5Q qualities.
        int res = readRecord(id, seq, quals, seqStream);
        if (res != 0)
            std::cerr << "Error reading file!\n";
     
        std::cout << id << '\t' << seq << '\t' << quals << '\n';
    }

    return 0;
}

