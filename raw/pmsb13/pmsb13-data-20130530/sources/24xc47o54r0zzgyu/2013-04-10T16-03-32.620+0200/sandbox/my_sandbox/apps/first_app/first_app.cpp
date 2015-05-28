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
	
	if (build(faiIndex, argv[1]) != 0)
	{
		std::cerr << "ERROR: Could not build the index!\n";
		return 1;
	}
	/*if (write(faiIndex, "C:/Users/Hox/Desktop/index.fasta.fai") != 0)
	{
		std::cerr << "ERROR: Could not write the index to file!\n";
		return 1;
	}*/
    return 0;
}