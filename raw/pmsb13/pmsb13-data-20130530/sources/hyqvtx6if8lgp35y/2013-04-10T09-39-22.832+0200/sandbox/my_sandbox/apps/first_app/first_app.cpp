#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align.h>

using namespace seqan;

typedef Align<DnaString,ArrayGaps> TAlign;
typedef Row<TAlign>::Type TRow;

int main()
{
	//Build sequences
	DnaString seq0 = "ACGTCACCTC";
	DnaString seq1 = "ACGGGCCTATC";
	//Build alignment data structure	
	TAlign align;
	resize(row(align), 2);
	assignSource(row(align,0),seq0);
    assignSource(row(align,1),seq1);
	//Insert gaps
	insertGaps(row(align,0),1,4);
	insertGaps(row(align,1),1,8);
	//Iterate over alignment to print the total count of gaps in it
	Iterator<, Rooted>::Type it = begin();
	unsigned count = 0;	
	for(; !atEnd(it); goNext(it))
	{
		if(isGap(it))
			++count;
	}
	std::cout << count << std::endl;
    return 0;
}
