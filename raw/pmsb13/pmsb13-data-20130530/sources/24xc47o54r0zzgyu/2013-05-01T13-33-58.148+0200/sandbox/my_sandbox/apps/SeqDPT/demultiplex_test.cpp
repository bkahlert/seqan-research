#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/find.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>

using namespace seqan;


SEQAN_DEFINE_TEST(findExactIndex_test)
{
	StringSet<String<Dna> > barcodes;
	append(barcodes, "AAAAAA");
	append(barcodes, "CCCCCC");
	append(barcodes, "GGGGGG");
	append(barcodes, "TTTTTT");
	append(barcodes, "ACGTAC");

	String<Dna5> readPiece = "CCCCCC"	

	Index<Dna5,IndexEsa<> > indices;
	Finder<esaIndex::Type, IndexEsa<> > finders;

	makeExactIndices(barcodes, indices, finders);
	res=findExactindex(readPiece, finders);
	SEQAN_ASSERT_EQ(1, res);
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
	SEQAN_CALL_TEST(findExactIndex_test); 
}
SEQAN_END_TESTSUITE
