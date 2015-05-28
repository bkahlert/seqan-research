#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include "demultiplex.h"

using namespace seqan;

SEQAN_DEFINE_TEST(check_test)
{
	StringSet<String<Dna5Q> > seqs;
	appendValue(seqs, "TACGTAGCTACTGACTGACT");
	appendValue(seqs, "G");
	appendValue(seqs, "TACGTCTGACT");
	appendValue(seqs, "");
	appendValue(seqs, "ATCGA");

	StringSet<String<char> > ids;
	appendValue(ids, "ErsteSeq");
	appendValue(ids, "LoeschenEins");
	appendValue(ids, "ZweiteSeq");
	appendValue(ids, "LoeschenLeer");
	appendValue(ids, "LoeschenFuenf");

	StringSet<String<Dna5Q> > barcodesFalse;
	appendValue(barcodesFalse, "ACGAGT");
	appendValue(barcodesFalse, "TGCATC");
	appendValue(barcodesFalse, "AGCTAAT");
	appendValue(barcodesFalse, "GTGACA");

	StringSet<String<Dna5Q> > barcodesTrue;
	appendValue(barcodesTrue, "ACGAGT");
	appendValue(barcodesTrue, "TGCATC");
	appendValue(barcodesTrue, "AGCTAA");
	appendValue(barcodesTrue, "GTGACA");

	StringSet<String<Dna5Q> > exspectedSeqs;
	appendValue(exspectedSeqs, "TACGTAGCTACTGACTGACT");
	appendValue(exspectedSeqs, "TACGTCTGACT");

	StringSet<String<char> > exspectedIds;
	appendValue(exspectedIds, "ErsteSeq");
	appendValue(exspectedIds, "ZweiteSeq");

	bool resF = check(seqs, ids, barcodesFalse);
	bool resT = check(seqs, ids, barcodesTrue);

	SEQAN_ASSERT_EQ(false, resF);
	SEQAN_ASSERT_EQ(true, resT);
	SEQAN_ASSERT_EQ(2, length(seqs));
	SEQAN_ASSERT_EQ(2, length(ids));
	SEQAN_ASSERT_EQ(exspectedSeqs[0], seqs[0]);
	SEQAN_ASSERT_EQ(exspectedIds[0], ids[0]);
	SEQAN_ASSERT_EQ(exspectedSeqs[1], seqs[1]);
	SEQAN_ASSERT_EQ(exspectedIds[1], ids[1]);
}

SEQAN_DEFINE_TEST(getPrefix_test) //Tests if the prefices are propperly extracted and too short sequences are excluded 
{
	StringSet<String<Dna5Q> > given;
	appendValue(given, "GATACAGACTGAGCATGTGATCGAC");
	//appendValue(given, "GAT"); //checked beforehand by different function (TODO)
	appendValue(given, "AATTCCGTACGTAGCTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	appendValue(given, "GTTGGAGTACGTAGCTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	
	StringSet<String<Dna5Q> > exspectedSet;
	appendValue(exspectedSet, "GATACAGACTGAGCATGTGATCGAC");
	appendValue(exspectedSet, "AATTCCGTACGTAGCTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	appendValue(exspectedSet, "GTTGGAGTACGTAGCTACGTACGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC");

	StringSet<String<Dna5Q> > exspected;
	appendValue(exspected, "GATACA");
	appendValue(exspected, "AATTCC");
	appendValue(exspected, "GTTGGA");
		
	StringSet<String<Dna5Q> > res = getPrefix(given, 6);
	for(unsigned i = 0; i < length(exspected); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], res[i]);
		SEQAN_ASSERT_EQ(exspectedSet[i], given[i]);
	}
}

SEQAN_DEFINE_TEST(findExactIndex_test) //Tests the exact search for a single prefix in the barcode index (therfore implicitly checks the construction of the index)
{
	StringSet<String<Dna5Q> > barcodes;
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "TTTTTT");
	appendValue(barcodes, "ACGTAC");

	StringSet<String<Dna5Q> > readPieces;
	appendValue(readPieces, "CCCCCC");
	appendValue(readPieces, "AAAAAA");
	appendValue(readPieces, "TTTTTT");
	appendValue(readPieces, "GGGGGG");
	appendValue(readPieces, "CCCNCC");
	appendValue(readPieces, "GATACA");
	appendValue(readPieces, "ACGTAC");
	appendValue(readPieces, "ATGACNAANG");	//can't happen in the first place...

	Index<StringSet<String<Dna5Q> >, IndexEsa<> > indexSet(barcodes);
	Finder<Index<StringSet<String<Dna5Q> >, IndexEsa<> > > esaFinder(indexSet);

	int exspected[] = {1,0,3,2,-1,-1,4,-1};

	for (unsigned i = 0; i < length(readPieces); ++i)
	{
		int res = findExactIndex(readPieces[i], esaFinder);
		SEQAN_ASSERT_EQ(exspected[i], res);
	}
}

SEQAN_DEFINE_TEST(findAllExactIndex_test) //Tests the exact search for a all prefices in the barcode index (therfore implicitly checks the construction of the index)
{
	StringSet<String<Dna5Q> > barcodes;
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "TTTTTT");
	appendValue(barcodes, "ACGTAC");
	
	StringSet<String<Dna5Q> > readPieces;
	appendValue(readPieces, "CCCCCC");
	appendValue(readPieces, "AAAAAA");
	appendValue(readPieces, "TTTTTT");
	appendValue(readPieces, "GGGGGG");
	appendValue(readPieces, "CCCNCC");
	appendValue(readPieces, "GATACA");
	appendValue(readPieces, "ACGTAC");

	Index<StringSet<String<Dna5Q> >, IndexEsa<> > indexSet(barcodes);
	Finder<Index<StringSet<String<Dna5Q> >, IndexEsa<> > > esaFinder(indexSet);

	int exspected[] = {1,0,3,2,-1,-1,4,};
	std::vector<int> res = findAllExactIndex(readPieces, esaFinder);
	for(unsigned i = 0; i < length(res); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], res[i]);
	}
	
}

SEQAN_DEFINE_TEST(findApprox_test) //Tests the approxmiate search for a single prefix in the barcodes
{
	StringSet<String<Dna5Q> > barcodes;
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "ACGTCA");
	appendValue(barcodes, "GATACA");
	std::vector<Pattern<String<Dna5Q>, DPSearch<SimpleScore> > > patterns = makePatterns(barcodes);
	
	StringSet<String<Dna5Q> > readPieces;
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
	StringSet<String<Dna5Q> > barcodes;
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "ACGTCA");
	appendValue(barcodes, "GATACA");
	std::vector<Pattern<String<Dna5Q>, DPSearch<SimpleScore> > > patterns = makePatterns(barcodes);

	StringSet<String<Dna5Q> > readPieces;
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
	StringSet<String<Dna5Q> > seqs;
	appendValue(seqs, "AAAAAAGTGACTGATCGTACGACTG");
	appendValue(seqs, "GGGGGGGGGGGGGGGG");

	std::vector<int> matches;
	appendValue(matches, 0);
	appendValue(matches, -1);

	StringSet<String<Dna5Q> > exspected ;
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
	StringSet<String<Dna5Q> > seqs;
	appendValue(seqs, "AAAAAAGTGACTGATCGTACGACTG");
	appendValue(seqs, "GGGGGGGGGGGGGGGG");

	StringSet<String<Dna5Q> > exspected ;
	appendValue(seqs, "GTGACTGATCGTACGACTG");
	appendValue(seqs, "GGGGGGGGGG");

	clipBarcodes(seqs, 6);
	for(unsigned i = 0; i < length(exspected); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], seqs[i]);
	}
	
}

SEQAN_DEFINE_TEST(group_test) //Tests if the grouping function returns correct groups. Implicitly tests resizeGroups
{
	StringSet<String<Dna5Q> > seqs;
	appendValue(seqs, "NANANANANANANANANANANANANANANANAN");
	appendValue(seqs, "NCNCNCNCNCNCNCNCNCNCNCNCN");
	appendValue(seqs, "NGNGNGNGNGNGNGNGNGNGNGNGNGNGNGNGN");
	appendValue(seqs, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
	appendValue(seqs, "GGGGGGGGGGGGG");
	appendValue(seqs, "NTNTNTNTNTNTNTNTNTNTNTNTNTNTNTNTNTNTNTNTNTNTNTN");

	StringSet<String<Dna5Q> > barcodes;
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "TTTTTT");
	appendValue(barcodes, "GAGAGA"); //barcode without matching sequences 
	appendValue(barcodes, "CCCCCC");
	
	std::vector<int> matches;
	appendValue(matches, 0);
	appendValue(matches, 4);
	appendValue(matches, 1);
	appendValue(matches, -1);
	appendValue(matches, 1);
	appendValue(matches, 2);

	std::vector<std::vector<String<Dna5Q>* > > exspected;
	resize(exspected, 7);
	appendValue(exspected[0], &(seqs[3]));
	appendValue(exspected[1], &(seqs[0]));
	appendValue(exspected[2], &(seqs[2]));
	appendValue(exspected[2], &(seqs[4]));	//expected[4] is left empty since the barcode has no matching sequences
	appendValue(exspected[3], &(seqs[5]));
	appendValue(exspected[5], &(seqs[1]));
	
	std::vector<std::vector<String<Dna5Q>* > > res = group(seqs, matches, barcodes);
	
	for(unsigned i = 0; i < length(exspected); ++i)
	{
		for(unsigned j = 0; j < length(exspected[i]); ++j)
		{
			SEQAN_ASSERT_EQ(exspected[i][j], res[i][j]);
		}
	}
}

SEQAN_DEFINE_TEST(DoAll_Exact_test)
{
	StringSet<String<Dna5Q> > seqs;
	appendValue(seqs, "AAAAAAGTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG");
	appendValue(seqs, "GGGGGGAGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");
	appendValue(seqs, "AAAAAATAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	appendValue(seqs, "GGGGGG");
	appendValue(seqs, "C");
	appendValue(seqs, "");
	appendValue(seqs, "TAGTCATACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");

	StringSet<String<char> > ids;
	appendValue(ids, "Adenin1");
	appendValue(ids, "Guanin");
	appendValue(ids, "Adenin2");
	appendValue(ids, "GuaninLoeschen");
	appendValue(ids, "CytosinLoeschen");
	appendValue(ids, "Loeschen");
	appendValue(ids, "Unidentifiziert");

	StringSet<String<char> > barcodes;
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "AAAAAA");

	StringSet<String<Dna5Q> > exspectedSeqs;
	appendValue(exspectedSeqs, "GTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG");
	appendValue(exspectedSeqs, "AGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");
	appendValue(exspectedSeqs, "TAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	appendValue(exspectedSeqs, "TAGTCATACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");

	std::vector<std::vector<String<Dna5Q>* > > exspected;
	resize(exspected, 4);
	appendValue(exspected[0], &(exspectedSeqs[3]));
	appendValue(exspected[1], &(exspectedSeqs[1])); //exspected[2] is left empty since the associated barcode is not matched
	appendValue(exspected[3], &(exspectedSeqs[0]));
	appendValue(exspected[3], &(exspectedSeqs[2]));

	std::vector<std::vector<String<Dna5Q>* > > res = DoAll(seqs, ids, barcodes, false, false);

	for(unsigned i = 0; i < length(exspected); ++i)
	{
		for(unsigned j = 0; j < length(exspected[i]); ++j)
		{
			SEQAN_ASSERT_EQ(*(exspected[i][j]), *(res[i][j]));
			std::cout<< barcodes[position(res[i][j])]<<std::endl;;
		}
	}
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
	SEQAN_CALL_TEST(check_test);
	SEQAN_CALL_TEST(getPrefix_test);
	SEQAN_CALL_TEST(findExactIndex_test); 
	SEQAN_CALL_TEST(findAllExactIndex_test); 
	SEQAN_CALL_TEST(findApprox_test);
	SEQAN_CALL_TEST(findAllApprox_test);
	SEQAN_CALL_TEST(clipBarcodes_test);
	SEQAN_CALL_TEST(clipBarcodesStrict_test);
	SEQAN_CALL_TEST(group_test);
	SEQAN_CALL_TEST(DoAll_Exact_test);
}
SEQAN_END_TESTSUITE