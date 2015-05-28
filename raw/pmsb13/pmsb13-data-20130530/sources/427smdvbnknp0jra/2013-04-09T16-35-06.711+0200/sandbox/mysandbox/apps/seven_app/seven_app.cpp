#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
    Dna5String genome = "ANTGGTTNCAACNGTAANTGCTGANNNACATGTNCGCGTGTA";
	Iterator<Dna5String >::Type it = begin(genome);
    for (goBegin(it,genome); !atEnd(it,genome); ++it){
        if (*it == 'N')
			*it = 'A';
    }
    std::cout << "Modified genome: " << genome << std::endl;
    return 0;
}