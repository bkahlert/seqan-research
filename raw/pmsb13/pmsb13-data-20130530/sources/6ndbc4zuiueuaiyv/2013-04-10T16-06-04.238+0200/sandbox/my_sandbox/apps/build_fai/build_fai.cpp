#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>


int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <fasta-file>" << std::endl;
        return 1;
    }
    
    seqan::FaiIndex faiIndex;
    int res = build(faiIndex, "/home/stefan/Dokumente/subversion/seqan-trunk/sandbox/my_sandbox/apps/build_fai/example.fasta");
    if (res != 0)
        std::cerr << "ERROR: Could not build the index!\n";

    std::cout << "faiIndex: " << faiIndex << std::endl;
    
    return 0;
}

