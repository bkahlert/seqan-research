#include <iostream>
#include <seqan/sequence.h>

using namespace seqan;

int main()
{
    Dna5String genome = "ANTGGTTNCAACNGTAANTGCTGANNNACATGTNCGCGTGTA";
    Iterator<Dna5String>::Type it = begin(genome);
    Iterator<Dna5String>::Type itEnd = end(genome);
    for (; it != itEnd; goNext(it)){
        if (genome[i] == 'N')
            value[i] = 'A';
    }
    std::cout << "Modified genome: " << genome << std::endl;
    return 0;
}