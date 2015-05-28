#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
    Dna5String genome = "ANTGGTTNCAACNGTAANTGCTGANNNACATGTNCGCGTGTA";
    for (Iterator<Dna5String,Rooted>::Type it1=begin(genome); !atEnd(it1); goNext(it1)){
	      if (getValue(it1) == (Dna5String)'N')
              assignValue(it1,'A');
    }
    std::cout << "Modified genome: " << genome << std::endl;
	system("Pause");
    return 0;
}