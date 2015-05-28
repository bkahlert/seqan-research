#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
    Dna5String genome = "ANTGGTTNCAACNGTAANTGCTGANNNACATGTNCGCGTGTA";
	Iterator <Dna5String, Standard>::Type it;

	for(it = begin(genome); it != end(genome); goNext(it)) {
		if((*it) == 'N') (*it) = 'A';
	}
	
    std::cout << "Modified genome: " << genome << std::endl;
    return 0;
}