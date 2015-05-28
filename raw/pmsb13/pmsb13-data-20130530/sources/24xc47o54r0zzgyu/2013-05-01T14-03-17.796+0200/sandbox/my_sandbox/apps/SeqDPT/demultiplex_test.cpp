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

	String<Dna> readPiece = "CCCCCC";	

	Index<StringSet<String<Dna> >, IndexEsa<> > indexSet(barcodes);
	Finder<Index<StringSet<String<Dna> >, IndexEsa<> > > esaFinder(indexSet);

	int res = findExactIndex(readPiece, esaFinder);
	SEQAN_ASSERT_EQ(1, res);
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
	SEQAN_CALL_TEST(findExactIndex_test); 
}
SEQAN_END_TESTSUITE
