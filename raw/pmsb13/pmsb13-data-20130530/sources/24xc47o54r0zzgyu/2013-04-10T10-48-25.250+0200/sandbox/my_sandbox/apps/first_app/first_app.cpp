#include <iostream>
#include <seqan/align.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>

using namespace seqan;




int main()
{
	Dna seq1 = "AAATGACGGATTG";
	Dna seq2 = "AGTCGGATCTACTG";

	typedef String<Dna> TSequence;
    typedef Align<TSequence,ArrayGaps> TAlign;

	TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0),seq1);
    assignSource(row(align,1),seq2);

	std::cout<<align;
	
	return 0;
}