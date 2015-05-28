#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
    Dna5String genome = "ANTGGTTNCAACNGTAANTGCTGANNNACATGTNCGCGTGTA";
	Iterator<Dna5String, Rooted >::Type it = begin(genome);
    for (; !atEnd(it); ++it){
        if (*it == 'N')
			*it = 'A';
    }
    std::cout << "Modified genome: " << genome << std::endl;
    return 0;
}