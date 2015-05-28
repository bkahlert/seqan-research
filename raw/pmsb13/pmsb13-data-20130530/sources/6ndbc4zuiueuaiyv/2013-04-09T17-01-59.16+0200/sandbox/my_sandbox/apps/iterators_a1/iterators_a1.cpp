#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
    Dna5String genome = "ANTGGTTNCAACNGTAANTGCTGANNNACATGTNCGCGTGTA";
    Iterator<Dna5String >::Type it = begin(genome);
    Iterator<Dna5String >::Type itEnd = end(genome);
    for (goBegin(it); !atEnd(it); goNext(it))
    {
        if (getValue(it) == 'N')
            genome(it) = 'A';
    }
    std::cout << "Modified genome: " << genome << std::endl;
    return 0;
}

