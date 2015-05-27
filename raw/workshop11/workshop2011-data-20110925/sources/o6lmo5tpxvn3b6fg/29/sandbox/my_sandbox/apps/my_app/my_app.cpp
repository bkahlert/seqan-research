#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout

typedef seqan::String<seqan::Dna> DnaString;

template <class A>
int count1mers(typename seqan::String<A> str)
{
	int alphabetsize = seqan::ValueSize<A>::VALUE;
	seqan::String<int> charcount;
	seqan::resize(charcount, alphabetsize, 0);

	for(int i = 0; )

	return 1;
}

int main(int, char **) {
	DnaString testseq = "CAGTGTGAGCGATTACATACCATCGACTAG";
    std::cout << testseq << std::endl;
	count1mers<seqan::Dna>(testseq);
    return 1;
}