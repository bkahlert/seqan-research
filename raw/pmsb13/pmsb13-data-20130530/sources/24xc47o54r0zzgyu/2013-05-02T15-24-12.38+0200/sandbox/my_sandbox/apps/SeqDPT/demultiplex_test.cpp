#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include "demultiplex.h"

using namespace seqan;

SEQAN_DEFINE_TEST(getPrefix_test)
{
	StringSet<String<Dna5> > given;
	appendValue(given, "GATACAGACTGAGCATGTGATCGAC");
	appendValue(given, "GAT");
	appendValue(given, "AATTCCGTACGTAGCTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	

	StringSet<String<Dna5> > exspected;
	appendValue(exspected, "GATACA");
	appendValue(exspected, "AATTCC");
		
	StringSet<String<Dna5> > res = getPrefix(given, 6);
	for(unsigned i = 0; i < length(exspected); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], res[i]);
	}
}

SEQAN_DEFINE_TEST(findExactIndex_test)
{
	StringSet<String<Dna> > barcodes;
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "TTTTTT");
	appendValue(barcodes, "ACGTAC");

	StringSet<String<Dna5> > readPieces;
	appendValue(readPieces, "CCCCCC");
	appendValue(readPieces, "AAAAAA");
	appendValue(readPieces, "TTTTTT");
	appendValue(readPieces, "GGGGGG");
	appendValue(readPieces, "CCCNCC");
	appendValue(readPieces, "GATACA");
	appendValue(readPieces, "ACGTAC");
	appendValue(readPieces, "ATGACNAANG"); //darf eigentlich eh nicht passieren

	Index<StringSet<String<Dna> >, IndexEsa<> > indexSet(barcodes);
	Finder<Index<StringSet<String<Dna> >, IndexEsa<> > > esaFinder(indexSet);

	int exspected[] = {1,0,3,2,-1,-1,4,-1};

	for (unsigned i = 0; i < length(readPieces); ++i)
	{
		int res = findExactIndex(readPieces[i], esaFinder);
		SEQAN_ASSERT_EQ(exspected[i], res);
	}
}

SEQAN_DEFINE_TEST(findAllExactIndex_test)
{
	StringSet<String<Dna> > barcodes;
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "TTTTTT");
	appendValue(barcodes, "ACGTAC");
	
	StringSet<String<Dna5> > readPieces;
	appendValue(readPieces, "CCCCCC");
	appendValue(readPieces, "AAAAAA");
	appendValue(readPieces, "TTTTTT");
	appendValue(readPieces, "GGGGGG");
	appendValue(readPieces, "CCCNCC");
	appendValue(readPieces, "GATACA");
	appendValue(readPieces, "ACGTAC");

	Index<StringSet<String<Dna> >, IndexEsa<> > indexSet(barcodes);
	Finder<Index<StringSet<String<Dna> >, IndexEsa<> > > esaFinder(indexSet);

	int exspected[] = {1,0,3,2,-1,-1,4,};
	std::vector<int> res = findAllExactIndex(readPieces, esaFinder);
	for(unsigned i = 0; i < length(res); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], res[i]);
	}
	
}

SEQAN_DEFINE_TEST(findApprox_test)
{
	String<Dna5> barcode = "CCCCCC";
	
	StringSet<String<Dna5> > readPieces;
	appendValue(readPieces, "CCCNCC");
	appendValue(readPieces, "CNCNCC");
	appendValue(readPieces, "CCCCCC");

	Finder<String<Dna5> > finder(barcode);
	
	int exspected[] = {1,0,0};

	for (unsigned i = 0; i < length(readPieces); ++i)
	{
		int res = findApprox(readPieces[i], barcode);
		SEQAN_ASSERT_EQ(exspected[i], res);
	}
}


SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
	SEQAN_CALL_TEST(findExactIndex_test); 
	SEQAN_CALL_TEST(findAllExactIndex_test); 
	SEQAN_CALL_TEST(getPrefix_test);
	SEQAN_CALL_TEST(findApprox_test);
}
SEQAN_END_TESTSUITE
