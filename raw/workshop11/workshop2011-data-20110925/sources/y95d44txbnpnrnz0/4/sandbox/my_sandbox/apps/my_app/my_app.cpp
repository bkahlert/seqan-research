#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <iostream>

using namespace seqan;


// Program entry point
int main()
	{
		typedef String<AminoAcid> TAminoAcidString; // define type aminoacid aplphabet plus string
		TAminoAcidString seq = "MQDRVKRPMNAFIVWSRDQRRKMALEN"; // declare the variable seq from type TAminoAcidString
		typedef Iterator<TAminoAcidString>::Type TIter; // define the type of iterator
		TIter itEnd = end(seq);
		
		typedef Size<TAminoAcidString>::Type TSize;
		typedef String<TSize> TCounterString;
		TCounterString counter;
		TSize alphSize = ValueSize<AminoAcid>::VALUE;
		resize(counter, alphSize, 0);
		
		for (TIter i = begin(seq); i != itEnd; goNext(i))
		{
			if (value(i) == 'R') value(i) = 'A';
			std::cout << value(i);
			value(counter, ordValue(value(i))) += 1;
		}
		std::cout << std::endl; //add a new line
		typedef Iterator<TCounterString>::Type TCounterIter;
		TCounterIter countIt = begin(counter);
		TCounterIter countItEnd = end(counter);
		for (TSize pos = 0; countIt != countItEnd; ++countIt, ++pos)
			std::cout << AminoAcid(pos) << ':' << value(countIt) << std::endl;
		
		
		std::cout << std::endl;
		return 0;
	}