#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include "demultiplex.h"

using namespace seqan;

int loadSeqs(char const * path, StringSet<String<char> >& ids, StringSet<String<Dna5Q> >& seqs)
{
	SequenceStream seqStream(path, SequenceStream::READ);
	unsigned records = 1000;									//number of records to be read
	resize(ids, records);
	resize(seqs, records);
	

	if (!isGood(seqStream))
	{
		std::cerr << "Error while opening the sequence-file.\n";
		return 1;
	}
	
	if(readBatch(ids, seqs, seqStream, records) != 0)		
	{
		std::cerr << "Error while reading the sequences.\n";
		return 1;
	}
	return 0;
}

int loadBarcodes(char const * path, StringSet<String<char> >& bcids, StringSet<String<Dna5Q> >& bcs)
{
	SequenceStream bcStream(path, SequenceStream::READ);
	
	if(!isGood (bcStream))
	{
		std::cerr << "Error while opening barcode-file.\n";
		return 1;
	}

	if(readAll(bcids, bcs, bcStream) != 0)
	{
		std::cerr << "Error while reading the barcodes.\n";
		return 1;
	}
	return 0;
}

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

	std::vector<std::vector<int> > exspected;
	resize(exspected, 7);
	appendValue(exspected[0], 3);
	appendValue(exspected[1], 0);
	appendValue(exspected[2], 2);
	appendValue(exspected[2], 4);	//expected[4] is left empty since the barcode has no matching sequences
	appendValue(exspected[3], 5);
	appendValue(exspected[5], 1);
	
	std::vector<std::vector<int> > res = group(seqs, matches, barcodes);
	
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

	std::vector<std::vector<int> > exspected;
	resize(exspected, 4);
	appendValue(exspected[0], 3);
	appendValue(exspected[1], 1);	//exspected[2] is left empty since the associated barcode is not matched
	appendValue(exspected[3], 0);
	appendValue(exspected[3], 2);

	std::vector<std::vector<int> > res = DoAll(seqs, ids, barcodes, false, false);

	for(unsigned i = 0; i < length(exspected); ++i)
	{
		for(unsigned j = 0; j < length(exspected[i]); ++j)
		{
			SEQAN_ASSERT_EQ(exspected[i][j], res[i][j]);
			SEQAN_ASSERT_EQ(exspectedSeqs[exspected[i][j]], seqs[res[i][j]]);
		}
	}
}

SEQAN_DEFINE_TEST(DoAll_Approx_test)
{
	StringSet<String<Dna5Q> > seqs;
	appendValue(seqs, "AAAAAAGTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG");
	appendValue(seqs, "GGNGNGAGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");
	appendValue(seqs, "NAAAAATAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	appendValue(seqs, "TAGTCATACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");

	StringSet<String<char> > ids;
	appendValue(ids, "Adenin0F");
	appendValue(ids, "Unidentifiziert2Fehler");
	appendValue(ids, "Adenin1F");
	appendValue(ids, "Unidentifiziert");

	StringSet<String<char> > barcodes;
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "AAAAAA");

	StringSet<String<Dna5Q> > exspectedSeqs;
	appendValue(exspectedSeqs, "GTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG");
	appendValue(exspectedSeqs, "GGNGNGAGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");
	appendValue(exspectedSeqs, "TAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	appendValue(exspectedSeqs, "TAGTCATACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");

	std::vector<std::vector<int> > exspected;
	resize(exspected, 4);
	appendValue(exspected[0], 1);
	appendValue(exspected[0], 3);	
	appendValue(exspected[3], 0);
	appendValue(exspected[3], 2);

	std::vector<std::vector<int> > res = DoAll(seqs, ids, barcodes, true, false);

	for(unsigned i = 0; i < length(exspected); ++i)
	{
		for(unsigned j = 0; j < length(exspected[i]); ++j)
		{
			SEQAN_ASSERT_EQ(exspected[i][j], res[i][j]);
			SEQAN_ASSERT_EQ(exspectedSeqs[exspected[i][j]], seqs[res[i][j]]);
		}
	}
}

SEQAN_DEFINE_TEST(DoAll_Approx_Multiplex_test)
{
	StringSet<String<Dna5Q> > seqs;
	appendValue(seqs, "GTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG");
	appendValue(seqs, "AGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");
	appendValue(seqs, "TAGCTAGCTAGCTAGCTAGCTAGCTAGC");

	StringSet<String<Dna5Q> > multiplex;
	appendValue(multiplex, "AAAAAA");
	appendValue(multiplex, "GNGNGG");
	appendValue(multiplex, "AACAAA");
	
	StringSet<String<char> > ids;
	appendValue(ids, "Adenin0F");
	appendValue(ids, "Unidentifiziert2Fehler");
	appendValue(ids, "Adenin1F");

	StringSet<String<Dna5Q> > barcodes;
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "AAAAAA");

	StringSet<String<Dna5Q> > exspectedSeqs;
	appendValue(exspectedSeqs, "GTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG");
	appendValue(exspectedSeqs, "AGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");
	appendValue(exspectedSeqs, "TAGCTAGCTAGCTAGCTAGCTAGCTAGC");

	std::vector<std::vector<int> > exspected;
	resize(exspected, 4);
	appendValue(exspected[0], 1);	
	appendValue(exspected[3], 0);
	appendValue(exspected[3], 2);

	std::vector<std::vector<int> > res = DoAll(seqs, ids, multiplex, barcodes, true);

	for(unsigned i = 0; i < length(exspected); ++i)
	{
		for(unsigned j = 0; j < length(exspected[i]); ++j)
		{
			SEQAN_ASSERT_EQ(exspected[i][j], res[i][j]);
			SEQAN_ASSERT_EQ(exspectedSeqs[exspected[i][j]], seqs[res[i][j]]);
		}
	}
}

SEQAN_DEFINE_TEST(DoAll_Exact_Multiplex_test)
{
	StringSet<String<Dna5Q> > seqs;
	appendValue(seqs, "GTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG");
	appendValue(seqs, "AGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");
	appendValue(seqs, "TAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	appendValue(seqs, "GGGG");

	StringSet<String<Dna5Q> > multiplex;
	appendValue(multiplex, "AAAAAA");
	appendValue(multiplex, "GGGGGG");
	appendValue(multiplex, "GGCCGG");
	appendValue(multiplex, "AAAAAA");

	StringSet<String<char> > ids;
	appendValue(ids, "Adenin");
	appendValue(ids, "Guanin");
	appendValue(ids, "GuaninUnidentifiziert");
	appendValue(ids, "GuaninKurz");

	StringSet<String<Dna5Q> > barcodes;
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "AAAAAA");

	StringSet<String<Dna5Q> > exspectedSeqs;
	appendValue(exspectedSeqs, "GTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG");
	appendValue(exspectedSeqs, "AGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");
	appendValue(exspectedSeqs, "TAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	appendValue(exspectedSeqs, "GGGG");

	std::vector<std::vector<int> > exspected;
	resize(exspected, 4);
	appendValue(exspected[0], 2);
	appendValue(exspected[1], 1);	//exspected[2] is left empty since the associated barcode is not matched
	appendValue(exspected[3], 0);
	appendValue(exspected[3], 3);

	std::vector<std::vector<int> > res = DoAll(seqs, ids, multiplex, barcodes, false);

	for(unsigned i = 0; i < length(exspected); ++i)
	{
		for(unsigned j = 0; j < length(exspected[i]); ++j)
		{
			SEQAN_ASSERT_EQ(exspected[i][j], res[i][j]);
			SEQAN_ASSERT_EQ(exspectedSeqs[exspected[i][j]], seqs[res[i][j]]);
		}
	}
}

SEQAN_DEFINE_TEST(Output_test)
{
	StringSet<String<Dna5Q> > seqs;
	appendValue(seqs, "AAAAAAGTACGATCGAAAAAAAAAAAAAAAAAAAAATGCTACGATGCTACG");
	appendValue(seqs, "GGGGGGAGTACGTACGGGGGGGGGGGGGAGCTAGCTAC");
	appendValue(seqs, "AAAAAATAGCTAGCTAGCTAAAAAAAAAAAAAAAGCTAGC");
	appendValue(seqs, "GGGGGG");
	appendValue(seqs, "C");
	appendValue(seqs, "");
	appendValue(seqs, "TAGTCATACGTACGTAGCNNNNNNNNNNNNNNNGCTAGCTAGCTAC");

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

	StringSet<String<char> > barcodeIds;
	appendValue(barcodeIds, "G-Barcode");
	appendValue(barcodeIds, "C-Barcode");
	appendValue(barcodeIds, "A-Barcode");

	StringSet<String<Dna5Q> > exspectedSeqs;
	appendValue(exspectedSeqs, "GTACGATCGAAAAAAAAAAAAAAAAAAAAATGCTACGATGCTACG");
	appendValue(exspectedSeqs, "AGTACGTACGGGGGGGGGGGGGAGCTAGCTAC");
	appendValue(exspectedSeqs, "TAGCTAGCTAGCTAAAAAAAAAAAAAAAGCTAGC");
	appendValue(exspectedSeqs, "TAGTCATACGTACGTAGCNNNNNNNNNNNNNNNGCTAGCTAGCTAC");

	std::vector<std::vector<int> > groups = DoAll(seqs, ids, barcodes, false, false);
	std::vector<StringSet<String<Dna5Q> > > gSeqs;
	std::vector<StringSet<String<char> > > gIds;
	buildSets(seqs, ids, groups, gSeqs, gIds);
	
	SEQAN_ASSERT_EQ(exspectedSeqs[0], gSeqs[2][0]);
	SEQAN_ASSERT_EQ(exspectedSeqs[1], gSeqs[1][0]);
	SEQAN_ASSERT_EQ(exspectedSeqs[2], gSeqs[2][1]);
	SEQAN_ASSERT_EQ(exspectedSeqs[3], gSeqs[0][0]);

	String<char> path =  "D:/seqan/out/";
	SEQAN_ASSERT_EQ(0, writeGroups(gSeqs, gIds, barcodeIds, groups, path));
}

SEQAN_DEFINE_TEST(Output_Paired_End_test)
{
	StringSet<String<Dna5Q> > seqs;
	appendValue(seqs, "AAAAAAGTACGATCGAAAAAAAAAAAAAAAAAAAAATGCTACGATGCTACG");
	appendValue(seqs, "GGGGGGAGTACGTACGGGGGGGGGGGGGAGCTAGCTAC");
	appendValue(seqs, "AAAAAATAGCTAGCTAGCTAAAAAAAAAAAAAAAGCTAGC");
	appendValue(seqs, "GGGGGG");
	appendValue(seqs, "C");
	appendValue(seqs, "");
	appendValue(seqs, "TAGTCATACGTACGTAGCNNNNNNNNNNNNNNNGCTAGCTAGCTAC");

	StringSet<String<Dna5Q> > seqsRev;
	appendValue(seqsRev, "NCTGATCGTACGTAAAAAAAAAAAAAAAAAAAAAACATGCTACG");
	appendValue(seqsRev, "NTCGTACGTAGCGGGGGGGGGGGGGGGGGATCGTACGTACGTACG");
	appendValue(seqsRev, "NTGCATGCTACTGAAAAAAAAAAAAAAAATCATGCTAGCTAGC");
	appendValue(seqsRev, "NAAAAAAAAAATGCATCGTAGCTAGCTAGCTAGC");
	appendValue(seqsRev, "C");
	appendValue(seqsRev, "");
	appendValue(seqsRev, "NGTACGTAGCTGNNNNNNNNNNNNNNNAGTCGATCGATCG");

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

	StringSet<String<char> > barcodeIds;
	appendValue(barcodeIds, "G-Barcode");
	appendValue(barcodeIds, "C-Barcode");
	appendValue(barcodeIds, "A-Barcode");

	StringSet<String<Dna5Q> > exspectedSeqs;
	appendValue(exspectedSeqs, "GTACGATCGAAAAAAAAAAAAAAAAAAAAATGCTACGATGCTACG");
	appendValue(exspectedSeqs, "AGTACGTACGGGGGGGGGGGGGAGCTAGCTAC");
	appendValue(exspectedSeqs, "TAGCTAGCTAGCTAAAAAAAAAAAAAAAGCTAGC");
	appendValue(exspectedSeqs, "TAGTCATACGTACGTAGCNNNNNNNNNNNNNNNGCTAGCTAGCTAC");

	StringSet<String<Dna5Q> > exspectedSeqsRev;
	appendValue(exspectedSeqsRev, "NCTGATCGTACGTAAAAAAAAAAAAAAAAAAAAAACATGCTACG");
	appendValue(exspectedSeqsRev, "NTCGTACGTAGCGGGGGGGGGGGGGGGGGATCGTACGTACGTACG");
	appendValue(exspectedSeqsRev, "NTGCATGCTACTGAAAAAAAAAAAAAAAATCATGCTAGCTAGC");
	appendValue(exspectedSeqsRev, "NGTACGTAGCTGNNNNNNNNNNNNNNNAGTCGATCGATCG");

	std::vector<std::vector<int> > groups = DoAll(seqs, seqsRev, ids, barcodes, false, false);
	std::vector<StringSet<String<Dna5Q> > > gSeqs;
	std::vector<StringSet<String<Dna5Q> > > gSeqsRev;
	std::vector<StringSet<String<char> > > gIds;
	buildSets(seqs, seqsRev, ids, groups, gSeqs, gSeqsRev, gIds);
	
	SEQAN_ASSERT_EQ(exspectedSeqs[0], gSeqs[2][0]);
	SEQAN_ASSERT_EQ(exspectedSeqs[1], gSeqs[1][0]);
	SEQAN_ASSERT_EQ(exspectedSeqs[2], gSeqs[2][1]);
	SEQAN_ASSERT_EQ(exspectedSeqs[3], gSeqs[0][0]);
	SEQAN_ASSERT_EQ(exspectedSeqsRev[0], gSeqsRev[2][0]);
	SEQAN_ASSERT_EQ(exspectedSeqsRev[1], gSeqsRev[1][0]);
	SEQAN_ASSERT_EQ(exspectedSeqsRev[2], gSeqsRev[2][1]);
	SEQAN_ASSERT_EQ(exspectedSeqsRev[3], gSeqsRev[0][0]);

	String<char> path =  "D:/seqan/out/";
	SEQAN_ASSERT_EQ(0, writeGroups(gSeqs, gSeqsRev, gIds, barcodeIds, groups, path));
}

SEQAN_DEFINE_TEST(IO_test)
{
	StringSet<String<char> > ids; 
	StringSet<String<Dna5Q> > seqs;
	StringSet<String<char> > bcids;
	StringSet<String<Dna5Q> > bcs;
	String<char> seqpath = "D:/seqan/test/seqs.fastq";
	String<char> bcpath = "D:/seqan/test/barcodes.fasta";

	int ls = loadSeqs(toCString(seqpath), ids, seqs);
	int bs = loadBarcodes(toCString(bcpath), bcids, bcs);

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
	SEQAN_CALL_TEST(DoAll_Approx_test);
	SEQAN_CALL_TEST(DoAll_Approx_Multiplex_test);
	SEQAN_CALL_TEST(DoAll_Exact_Multiplex_test);
	SEQAN_CALL_TEST(Output_test);
	SEQAN_CALL_TEST(Output_Paired_End_test);
	SEQAN_CALL_TEST(IO_test);
}
SEQAN_END_TESTSUITE