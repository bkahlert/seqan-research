#include <iostream>
#include <vector>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout

typedef seqan::String<seqan::Dna> DnaString;

template <typename A>
int count1mers(seqan::String<A> str)
{
	std::vector<int> charcount (ValueSize<A>::VALUE,0);

	for(int i = 0; i < length(str); i++)
	{
		charcount[str[i]] ++;
	}

	for(int i = 0; i<256; i++)
	{
		if(charcount[i] != 0)
			std::cout << char(i) << " " << charcount[i] << std::endl;
	}
	return 1;
}

int main(int, char **) {
	DnaString testseq = "CAGTGTGAGCGATTACATACCATCGACTAGT";
    std::cout << testseq << std::endl;
	count1mers<seqan::Dna>(testseq);
    return 1;
}