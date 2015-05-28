#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/basic.h>

using namespace seqan;

int main()
{
    Dna5String genome = "ANTGGTTNCAACNGTAANTGCTGANNNACATGTNCGCGTGTA";
	Iterator <Dna5String, Standard>::Type it;

	for(it = begin(genome); it != end(genome); goNext(it)) {
		if((*it) == 'N') (*it) = 'A';
	}
	
    std::cout << "Modified genome: " << genome << std::endl;

	Iterator <Dna5String, Rooted>::Type it2 = begin(genome);
	//if(container(it2))

	for(it2 = begin(genome); it2 != end(genome); goNext(it2)) {
		if((*it2) == 'N') (*it2) = 'A';
	}
	
    std::cout << "Modified genome: " << genome << std::endl;

    return 0;
}