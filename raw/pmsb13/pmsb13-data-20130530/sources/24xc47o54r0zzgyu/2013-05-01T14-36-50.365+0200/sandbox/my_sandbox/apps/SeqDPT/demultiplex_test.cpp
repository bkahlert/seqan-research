#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/find.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include "demultiplex.h"

using namespace seqan;


SEQAN_DEFINE_TEST(findExactIndex_test)
{
	StringSet<String<Dna> > barcodes;
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "TTTTTT");
	appendValue(barcodes, "ACGTAC");

	StringSet<Dna> readPieces;
	appendValue(readPieces, "CCCCCC");
	appendValue(readPieces, "AAAAAA");
	appendValue(readPieces, "TTTTTT");
	appendValue(readPieces, "GGGGGG");
	appendValue(readPieces, "CCCCCC");
	appendValue(readPieces, "GATACA");
	appendValue(readPieces, "ACGTAC");

	Index<StringSet<String<Dna> >, IndexEsa<> > indexSet(barcodes);
	Finder<Index<StringSet<String<Dna> >, IndexEsa<> > > esaFinder(indexSet);

	std::vector<int> exspected = {1,0,3,2,1,-1,4};

	for (unsigned i = 0; i < length(readPieces); ++i)
	{
		int res = findExactIndex(readPieces[i], esaFinder);
		SEQAN_ASSERT_EQ(exspected, res);
	}
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
	SEQAN_CALL_TEST(findExactIndex_test); 
}
SEQAN_END_TESTSUITE
