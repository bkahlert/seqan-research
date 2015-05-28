#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/basic.h>
#include <seqan/file.h>
using namespace seqan;
int main(int argc, char const ** argv)
{
    if (argc < 2)
    {
        std::cerr << "USAGE: basic_seq_io_example FILENAME\n";
        return 1;
    }
    seqan::FaiIndex faiIndex;
    int res = build(faiIndex, argv[1]);
    if (res != 0)
        std::cerr << "ERROR: Could not build the index!\n";
    CharString tmp = argv[1];
    Append(tmp, ".fai");
    res = write(faiIndex, tmp);
    if (res != 0)
        std::cerr << "ERROR: Could not write the index to file!\n";

}

