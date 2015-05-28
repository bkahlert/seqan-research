#include <iostream>
#include <seqan/seq_io.h>

int main(int argc, char const ** argv)
{
    if (argc < 2)
    {
        std::cerr << "USAGE: basic_seq_io_example FILENAME\n";
        return 1;
    }

   
	seqan::FaiIndex faiIndex;
	int res = build(faiIndex, argv[1], "C:/Users/Hox/Desktop/index.fasta.fai");
	if (res != 0)
		std::cerr << "ERROR: Could not build the index!\n";
    return 0;
}