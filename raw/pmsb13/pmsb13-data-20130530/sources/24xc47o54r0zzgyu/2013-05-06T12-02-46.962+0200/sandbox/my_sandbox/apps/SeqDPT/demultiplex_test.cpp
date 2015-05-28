#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include "demultiplex.h"

using namespace seqan;

SEQAN_DEFINE_TEST(getPrefix_test) //Tests if the prefices are propperly extracted and too short sequences are excluded 
{
	StringSet<String<Dna5> > given;
	appendValue(given, "GATACAGACTGAGCATGTGATCGAC");
	//appendValue(given, "GAT"); //checked beforehand by different function (TODO)
	appendValue(given, "AATTCCGTACGTAGCTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	appendValue(given, "GTTGGAGTACGTAGCTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	
	StringSet<String<Dna5> > exspectedSet;
	appendValue(exspectedSet, "GATACAGACTGAGCATGTGATCGAC");
	appendValue(exspectedSet, "AATTCCGTACGTAGCTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	appendValue(exspectedSet, "GTTGGAGTACGTAGCTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC");

	StringSet<String<Dna5> > exspected;
	appendValue(exspected, "GATACA");
	appendValue(exspected, "AATTCC");
	appendValue(exspected, "GTTGGA");
		
	StringSet<String<Dna5> > res = getPrefix(given, 6);
	for(unsigned i = 0; i < length(exspected); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], res[i]);
		SEQAN_ASSERT_EQ(exspectedSet[i], given[i]);
	}
}

SEQAN_DEFINE_TEST(findExactIndex_test) //Tests the exact search for a single prefix in the barcode index (therfore implicitly checks the construction of the index)
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
	appendValue(readPieces, "ATGACNAANG");	//can't happen in the first place...

	Index<StringSet<String<Dna> >, IndexEsa<> > indexSet(barcodes);
	Finder<Index<StringSet<String<Dna> >, IndexEsa<> > > esaFinder(indexSet);

	int exspected[] = {1,0,3,2,-1,-1,4,-1};

	for (unsigned i = 0; i < length(readPieces); ++i)
	{
		int res = findExactIndex(readPieces[i], esaFinder);
		SEQAN_ASSERT_EQ(exspected[i], res);
	}
}

SEQAN_DEFINE_TEST(findAllExactIndex_test) //Tests the exact search for a all prefices in the barcode index (therfore implicitly checks the construction of the index)
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

SEQAN_DEFINE_TEST(findApprox_test) //Tests the approxmiate search for a single prefix in the barcodes
{
	StringSet<String<Dna5> > barcodes;
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "ACGTCA");
	appendValue(barcodes, "GATACA");
	std::vector<Pattern<String<Dna>, DPSearch<SimpleScore> > > patterns = makePatterns(barcodes);
	
	StringSet<String<Dna5> > readPieces;
	appendValue(readPieces, "CCCNCC");
	appendValue(readPieces, "CNCNCC");
	appendValue(readPieces, "CCCCC");
	appendValue(readPieces, "CCCNCC");
	appendValue(readPieces, "AAAAAA");
	appendValue(readPieces, "NAAAAA");
	appendValue(readPieces, "AAAAA");

	int exspected[] = {0, -1, -1, 0, 1, 1,-1};

	for(unsigned i = 0; i < length(readPieces); ++i)
	{
		int res = findApprox(readPieces[i], patterns);
		SEQAN_ASSERT_EQ(exspected[i], res);
	}
	
	
}

SEQAN_DEFINE_TEST(findAllApprox_test) //Tests the approxmiate search for all prefices in the barcodes
{
	StringSet<String<Dna5> > barcodes;
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "ACGTCA");
	appendValue(barcodes, "GATACA");
	std::vector<Pattern<String<Dna>, DPSearch<SimpleScore> > > patterns = makePatterns(barcodes);

	StringSet<String<Dna5> > readPieces;
	appendValue(readPieces, "CCCNCC");
	appendValue(readPieces, "CNCNCC");
	appendValue(readPieces, "CCCCC");
	appendValue(readPieces, "CCCNCC");
	appendValue(readPieces, "AAAAAA");
	appendValue(readPieces, "NAAAAA");
	appendValue(readPieces, "AAAAA");

	int exspected[] = {0, -1, -1, 0, 1, 1,-1};
	
	std::vector<int> res = findAllApprox(readPieces, patterns);
	for(unsigned i = 0; i < length(res); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], res[i]);
	}
} 

SEQAN_DEFINE_TEST(clipBarcodes_test) //Tests if barcodes are propperly clipped
{
	StringSet<String<Dna5> > seqs;
	appendValue(seqs, "AAAAAAGTGACTGATCGTACGACTG");
	appendValue(seqs, "GGGGGGGGGGGGGGGG");

	std::vector<int> matches;
	appendValue(matches, 0);
	appendValue(matches, -1);

	StringSet<String<Dna5> > exspected ;
	appendValue(seqs, "GTGACTGATCGTACGACTG");
	appendValue(seqs, "GGGGGGGGGGGGGGGG");

	clipBarcodes(seqs, matches, 6);
	for(unsigned i = 0; i < length(exspected); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], seqs[i]);
	}
	
}

SEQAN_DEFINE_TEST(clipBarcodesStrict_test) //Tests if the barcodes are propperly clipped using the strict clipping overload of clipBarcode
{
	StringSet<String<Dna5> > seqs;
	appendValue(seqs, "AAAAAAGTGACTGATCGTACGACTG");
	appendValue(seqs, "GGGGGGGGGGGGGGGG");

	StringSet<String<Dna5> > exspected ;
	appendValue(seqs, "GTGACTGATCGTACGACTG");
	appendValue(seqs, "GGGGGGGGGG");

	clipBarcodes(seqs, 6);
	for(unsigned i = 0; i < length(exspected); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], seqs[i]);
	}
	
}

SEQAN_DEFINE_TEST(group_test) //Tests if the grouping function return correct groups TODO: Fehlerhaft
{
	StringSet<String<Dna5> > seqs;
	appendValue(seqs, "NANANANANANANANANANANANANANANANAN");
	appendValue(seqs, "NCNCNCNCNCNCNCNCNCNCNCNCN");
	appendValue(seqs, "NGNGNGNGNGNGNGNGNGNGNGNGNGNGNGNGN");
	appendValue(seqs, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
	appendValue(seqs, "GGGGGGGGGGGGG");
	appendValue(seqs, "NTNTNTNTNTNTNTNTNTNTNTNTNTNTNTNTNTNTNTNTNTNTNTN");

	StringSet<String<char> > ids;
	appendValue(ids, "Adenin");
	appendValue(ids, "Cytosin");
	appendValue(ids, "Guanin");
	appendValue(ids, "Muell");
	appendValue(ids, "Mehr_Guanin");
	appendValue(ids, "Thymin");

	StringSet<String<Dna> > barcodes;
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "TTTTTT");
	appendValue(barcodes, "CCCCCC");
	
	std::vector<int> matches;
	appendValue(matches, 0);
	appendValue(matches, 3);
	appendValue(matches, 1);
	appendValue(matches, -1);
	appendValue(matches, 1);
	appendValue(matches, 2);

	seqAtts ahs;
	ahs.seq = seqs[0];
	ahs.id = ids[0];
	ahs.barcode = barcodes[0];

	seqAtts cehs;
	cehs.seq = seqs[1];
	cehs.id = ids[1];
	cehs.barcode = barcodes[3];

	seqAtts gehs;
	gehs.seq = seqs[2];
	gehs.id = ids[2];
	gehs.barcode = barcodes[1];

	seqAtts muell;
	muell.seq = seqs[3];
	muell.id = ids[3];
	muell.barcode = "NONE";

	seqAtts mehrgehs;
	mehrgehs.seq = seqs[4];
	mehrgehs.id = ids[4];
	mehrgehs.barcode = barcodes[1];

	seqAtts tehs;
	tehs.seq = seqs[5];
	tehs.id = ids[5];
	tehs.barcode = barcodes[2];
	
	std::vector<std::vector<seqAtts> > exspected;
	resize(exspected, 6);
	appendValue(exspected[0], ahs);
	appendValue(exspected[1], gehs);
	appendValue(exspected[1], mehrgehs);
	appendValue(exspected[2], tehs);
	appendValue(exspected[3], cehs);
	
	std::vector<std::vector<seqAtts> > res = group(seqs, ids, matches, barcodes);
	
	for(unsigned i = 0; i < length(exspected); ++i)
	{
		for(unsigned j = 0; j < length(exspected[i]); ++j)
		{
			SEQAN_ASSERT_EQ(exspected[i][j].seq, res[i][j].seq);
			SEQAN_ASSERT_EQ(exspected[i][j].id, res[i][j].id);
			SEQAN_ASSERT_EQ(exspected[i][j].barcode, res[i][j].barcode);
		}
	}
}


SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
	SEQAN_CALL_TEST(getPrefix_test);
	SEQAN_CALL_TEST(findExactIndex_test); 
	SEQAN_CALL_TEST(findAllExactIndex_test); 
	SEQAN_CALL_TEST(findApprox_test);
	SEQAN_CALL_TEST(findAllApprox_test);
	SEQAN_CALL_TEST(clipBarcodes_test);
	SEQAN_CALL_TEST(clipBarcodesStrict_test);
	SEQAN_CALL_TEST(group_test);
}
SEQAN_END_TESTSUITE