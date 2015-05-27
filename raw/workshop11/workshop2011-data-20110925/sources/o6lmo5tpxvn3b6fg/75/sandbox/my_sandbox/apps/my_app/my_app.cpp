#include <iostream>
#include <vector>

#include <seqan/align.h>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout

typedef seqan::String<seqan::Dna> DnaString;
typedef seqan::String<seqan::AminoAcid> PeptideString;
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
			std::cout << A(i) << " " << charcount[i] << std::endl;
	}
	return 1;
}

void allStrings_append(CharString s, int restlen)
{
	if (restlen == 0)
	{
		std::cout << s << std::endl;
	} 
	else
	{
		for(char c = 'a'; c <= 'z'; c++)
		{
			CharString temp = s;
			seqan::appendValue(temp,c);
			allStrings_append(temp,restlen-1);
		}
	}
}


int allStrings(int len)
{
	allStrings_append("",len);
	return 1;
}


void replaceX(CharString str)
{
	for(auto it = seqan::begin(str); it!=seqan::end(str); it++)
	{
		if(*it=='a')
			(*it) = 'X';
	}
	std::cout << str << std::endl;

}


void align(DnaString s1, DnaString s2)
{
	typedef seqan::Align<DnaString, seqan::ArrayGaps> TAlign;
	TAlign align;
	seqan::resize(seqan::rows(align),2);
	seqan::assignSource(seqan::row(align,0),s1);
	seqan::assignSource(seqan::row(align,1),s2);
	int score = seqan::globalAlignment(align,seqan::Score<int>(1,-1,-1,-1), seqan::Hirschberg());

	std::cout << "Alignmnent score: " << score << std::endl;

	auto it1 = seqan::begin(seqan::row(align,0));
	auto end1 = seqan::end(seqan::row(align,0));
	auto it2 = seqan::begin(seqan::row(align,1));
	auto end2 = seqan::end(seqan::row(align,1));

	for(; it1 != end1; it1++)
		if(seqan::isGap(it1))
			std::cout << seqan::position(it1) << ", ";
	std::cout << std::endl;

	
}

int main(int, char **) {
	//DnaString testseq = "CAGTGTGAGCGATTACATACCATCGACTAGT";
 //   std::cout << testseq << std::endl;
	//count1mers(testseq);

	//PeptideString testpep = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
 //   std::cout << testpep << std::endl;
	//count1mers(testpep);

	//CharString testchar = "Hello World!";
 //   std::cout << testchar << std::endl;
	//count1mers(testchar);

//	allStrings(3);

	replaceX("bla");



    return 1;
}