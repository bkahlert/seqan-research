#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align.h>

using namespace seqan;

//Define all types
typedef DnaString TSequence;
typedef Align<TSequence,ArrayGaps> TAlign;
typedef Row<TAlign>::Type TRow;
typedef Iterator<TRow>::Type TRowIterator;

int main()
{
	//Build sequences
	TSequence seq0 = "ACGTCACCTC";
	TSequence seq1 = "ACGGGCCTATC";
	//Build alignment data structure	
	TAlign align;
	resize(rows(align), 2);
	assignSource(row(align,0),seq0);
    assignSource(row(align,1),seq1);
	//Insert gaps
	TRow & row0 = row(align,0);
    TRow & row1 = row(align,1);
	insertGaps(row0,2,2);
    insertGap(row0,7);
    insertGaps(row1,9,2);
	//Iterate over alignment to print the total count of gaps in it
	TRowIterator itRow0 = begin(row0);
    TRowIterator itEndRow0 = end(row0);
    TRowIterator itRow1 = begin(row1);
	int count = 0;	
	for(; !atEnd(itRow0, row0); goNext(itRow0), goNext(itRow1))
    {
        if(isGap(itRow0))
        {
            ++count;
        }
        if(isGap(itRow1))
        {
            ++count;
        }
    }
	std::cout << "Gaps: " << count << std::endl;
	std::cout << align << std::endl;
    return 0;

