#include <iostream>
#include <seqan/align.h>

using namespace seqan;

//Define all types
typedef String<char> TSequence;
typedef Align<TSequence, ArrayGaps> TAlign;
typedef StringSet<TSequence> TStringSet;
typedef Gaps<TSequence, ArrayGaps> TGaps;
typedef Iterator<TGaps>::Type TGapsIterator;
typedef Iterator<String<int> >::Type TIterator;

int main()
{
	//Build strings	
	TSequence text = "MISSISSIPPIANDMISSOURI";
	TSequence pattern = "SISSI"; 
	//Build alignment data structure	
	TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0), text);
    assignSource(row(align,1), pattern);
	//Report the positions of all hits with at most 2 errors
	String<int> positions;
	for(unsigned i = 0; i < length(text)-length(pattern); ++i)
	{
		TSequence new_text = infix(text, i, i+length(pattern));
		if(globalAlignmentScore(new_text, pattern, MyersHirschberg())<2)
			appendValue(positions, i);	
	}
	std::cout << positions << std::endl;
	return 0;
}
