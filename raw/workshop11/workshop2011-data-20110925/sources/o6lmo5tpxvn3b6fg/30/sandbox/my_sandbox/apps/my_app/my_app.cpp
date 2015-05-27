#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout

typedef seqan::String<seqan::Dna> DnaString;

int count1mers(seqan::String<char> str)
{
	seqan::String<int> charcount;
	seqan::resize(charcount, 256, 0);

	for(int i = 0; sizeof(str); i++)
	{
		charcount[char(str[i])] ++;
	}

	for(int i = 0; sizeof(charcount); i++)
	{
		if(charcount[i] != 0)
			std::cout << char(i) << " " << charcount[i] << std::endl;
	}
	return 1;
}

int main(int, char **) {
	DnaString testseq = "CAGTGTGAGCGATTACATACCATCGACTAG";
    std::cout << testseq << std::endl;
	count1mers<seqan::Dna>(testseq);
    return 1;
}