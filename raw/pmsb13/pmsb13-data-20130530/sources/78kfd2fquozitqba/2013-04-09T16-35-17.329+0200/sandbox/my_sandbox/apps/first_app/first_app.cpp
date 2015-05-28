#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
    Dna5String genome = "ANTGGTTNCAACNGTAANTGCTGANNNACATGTNCGCGTGTA";
    Iterator<Dna5String>::Type it = begin(genome);

    for (;it!=atEnd(genome);goNext(it)){
        if (genome[it] == 'N')
            genome[it] = 'A';
    }
    std::cout << "Modified genome: " << genome << std::endl;
    return 0;
}
