#include <iostream>
#include <vector>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout

typedef seqan::String<seqan::Dna> DnaString;
typedef seqan::String<seqan::Peptide> PeptideString;
typedef seqan::String<char> CharString;

template <typename A>
int count1mers(seqan::String<A> str)
{
	std::vector<int> charcount (seqan::ValueSize<A>::VALUE,0);

	for(int i = 0; i < length(str); i++)
	{
		charcount[seqan::ordValue(str[i])] ++;
	}

	for(int i = 0; i<seqan::ValueSize<A>::VALUE; i++)
	{
		if(charcount[i] != 0)
			std::cout << seqan::Dna(i) << " " << charcount[i] << std::endl;
	}
	return 1;
}

int main(int, char **) {
	DnaString testseq = "CAGTGTGAGCGATTACATACCATCGACTAGT";
    std::cout << testseq << std::endl;
	count1mers(testseq);

	PeptideString testpep = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
    std::cout << testpep << std::endl;
	count1mers(testpep);

	CharString testchar = "Hello World!";
    std::cout << testchar << std::endl;
	count1mers(testchar);
    return 1;
}