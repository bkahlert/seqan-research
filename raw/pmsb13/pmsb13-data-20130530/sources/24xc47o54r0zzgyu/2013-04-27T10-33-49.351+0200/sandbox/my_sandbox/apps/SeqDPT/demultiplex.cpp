/*! Program part for barcode demultiplexing */

#include <seqan/find.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>

using namespace seqan;

typedef Dna5String TSequence;
Index<TSequence> index;
Finder<Index<CharString> > finder;
StringSet<TSequence> needles;

/*!
 Function for concatenating the barcordes.
*/
String<TSequence> concatBc(const SequenceStream & bcs)
{

}


/*!
 Function for creating the haystacks of barcodes.
*/
int makeExactIndex(const TSequence & barcodes)
{
	return 0;
}

/*!
 Function for searching 1 piece of sequence in the barcode index.
*/
int findExactIndex(const TSequence & seq, const Index<TSequence> & bchaystack)
{
	return 0;
}

int main()
{
	return 0;
}