#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout

typedef seqan::String<seqan::Dna> DnaString;

template <A>
int count1mers(seqan::String<A> str)
{
	
	seqan::String<int> charcount;
	seqan::resize(charcount, 26, 0);


	return 1;
}

int main(int, char **) {
	DnaString testseq = "CAGTGTGAGCGATTACATACCATCGACTAG";
    std::cout << testseq << std::endl;

    return 1;
}