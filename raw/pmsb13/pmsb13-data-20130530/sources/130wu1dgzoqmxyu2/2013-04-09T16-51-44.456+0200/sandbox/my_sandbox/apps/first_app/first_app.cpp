#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
    Dna5String genome = "ANTGGTTNCAACNGTAANTGCTGANNNACATGTNCGCGTGTA";
    typedef Iterator<Dna5String>::Type TIterator;
    for (TIterator it = begin(genome); !atEnd(it, genome); ++i){
        if (value(it) == 'N')
        	value(it) = 'A';
    }
    std::cout << "Modified genome: " << genome << std::endl;
    return 0;
}
