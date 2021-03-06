// FRAGMENT(create-string)
#include <iostream>
#include <seqan/basic.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
	typedef String<AminoAcid> TAminoAcidString;
	TAminoAcidString sourceSeq = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
// FRAGMENT(iterate-and-replace)
	typedef Iterator<TAminoAcidString>::Type TIter;

	TIter itEnd = end(sourceSeq);
	for (TIter it = begin(sourceSeq); it != itEnd; goNext(it))
	{
		if (value(it) == 'R') value(it) = 'A';
		std::cout << value(it) << ',';
	}
	std::cout << std::endl;
// FRAGMENT(count-occurences)
	typedef Size<TAminoAcidString>::Type TSize;
	typedef String<TSize> TCounterString;
	TCounterString counter;
	TSize alphSize = ValueSize<AminoAcid>::VALUE;
	resize(counter, alphSize, 0);
	for (TIter it = begin(sourceSeq); it != itEnd; goNext(it))
		value(counter, ordValue(value(it))) += 1;
// FRAGMENT(frequency-table)
	typedef Iterator<TCounterString>::Type TCounterIter;
	TCounterIter countIt = begin(counter);
	TCounterIter countItEnd = end(counter);
	for (TSize pos = 0; countIt != countItEnd; ++countIt, ++pos)
		std::cout << AminoAcid(pos) << ':' << value(countIt) << std::endl;

	return 0;
}

