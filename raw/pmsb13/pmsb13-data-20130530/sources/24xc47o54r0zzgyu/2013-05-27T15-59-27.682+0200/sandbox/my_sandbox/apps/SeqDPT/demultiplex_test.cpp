#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include "demultiplex.h"

using namespace seqan;

// Used for loading the sequences. Needed for the io-test.
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
// Used for loading the sequences. Needed for the io-test.
int loadBarcodes(char const * path, StringSet<String<char> >& bcids, StringSet<String<Dna> >& bcs)
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
// Checks the correctness of the check function which checks the size of the barcodes and reads.
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
	SEQAN_ASSERT_EQ(2u, length(seqs));
	SEQAN_ASSERT_EQ(2u, length(ids));
	SEQAN_ASSERT_EQ(exspectedSeqs[0], seqs[0]);
	SEQAN_ASSERT_EQ(exspectedIds[0], ids[0]);
	SEQAN_ASSERT_EQ(exspectedSeqs[1], seqs[1]);
	SEQAN_ASSERT_EQ(exspectedIds[1], ids[1]);
}
// Checks the correctness of the getPrefix function which extract the prefices from a set of sequences.
SEQAN_DEFINE_TEST(getPrefix_test)
{
	StringSet<String<Dna5Q> > given;
	appendValue(given, "GATACAGACTGAGCATGTGATCGAC");
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
// Checks the correctness of the findExactIndex function which searches for one piece of sequence in the barcodes. Implicitly checks the construction of the Index.
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
// Checks the correctnes of the findAllExactIndex function which searches for many pieces of sequence in the barcodes. Implicitly checks the construction of the Index.
SEQAN_DEFINE_TEST(findAllExactIndex_test) 
{
	StringSet<String<Dna5Q> > barcodes;
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "TTTTTT");
	appendValue(barcodes, "ACGTAC");
	
	DemultiplexStats demultiplexStats;
	resize(demultiplexStats.groups, length(barcodes)+1);

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
	std::vector<int> res = findAllExactIndex(readPieces, esaFinder, demultiplexStats);
	for(unsigned i = 0; i < length(res); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], res[i]);
	}
	
}
// Checks the correctnes of the findApprox function which searches for one piece of sequence in the barcodes allowing one mismatch.
SEQAN_DEFINE_TEST(findApprox_test)
{
	StringSet<String<Dna5Q> > barcodes;
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "ACGTCA");
	appendValue(barcodes, "GATACA");
	std::vector<Pattern<String<Dna5Q>, DPSearch<SimpleScore> > > patterns = makePatterns(barcodes);

	DemultiplexStats demultiplexStats;
	resize(demultiplexStats.groups, length(barcodes)+1);
	
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
// Checks the correctnes of the findAllApprox function which searches for many pieces of sequence in the barcodes allowing one mismatch. 
SEQAN_DEFINE_TEST(findAllApprox_test)
{
	StringSet<String<Dna5Q> > barcodes;
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "AAAAAA");
	appendValue(barcodes, "ACGTCA");
	appendValue(barcodes, "GATACA");
	std::vector<Pattern<String<Dna5Q>, DPSearch<SimpleScore> > > patterns = makePatterns(barcodes);

	DemultiplexStats demultiplexStats;
	resize(demultiplexStats.groups, length(barcodes)+1);

	StringSet<String<Dna5Q> > readPieces;
	appendValue(readPieces, "CCCNCC");
	appendValue(readPieces, "CNCNCC");
	appendValue(readPieces, "CCCCC");
	appendValue(readPieces, "CCCNCC");
	appendValue(readPieces, "AAAAAA");
	appendValue(readPieces, "NAAAAA");
	appendValue(readPieces, "AAAAA");

	int exspected[] = {0, -1, -1, 0, 1, 1,-1};
	
	std::vector<int> res = findAllApprox(readPieces, patterns, demultiplexStats);
	for(unsigned i = 0; i < length(res); ++i)
	{
		SEQAN_ASSERT_EQ(exspected[i], res[i]);
	}
} 
// Checks the correctnes of the clipBarcodes function which erases the first x bases of a sequence.
SEQAN_DEFINE_TEST(clipBarcodes_test)
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
// Checks the correctnes of the clipBarcodes function using the hardClip method and therefore also clipping sequences without matching barcode.
SEQAN_DEFINE_TEST(clipBarcodesStrict_test)
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
// Checks the correctnes of the group function which builds the groups of barcodes matching the same barcodes. Implicitly checks resizeGroups.
SEQAN_DEFINE_TEST(group_test)
{
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
	
	std::vector<std::vector<int> > res = group(matches, barcodes);
	
	for(unsigned i = 0; i < length(exspected); ++i)
	{
		for(unsigned j = 0; j < length(exspected[i]); ++j)
		{
			SEQAN_ASSERT_EQ(exspected[i][j], res[i][j]);
		}
	}
}
// Checks the correctness of the DoAll function performing all demultiplexing operations for exact inline barcode matching.
SEQAN_DEFINE_TEST(DoAll_Exact_test)
{
	StringSet<String<Dna5Q> > seqs;
	appendValue(seqs, "AAAAAAGTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG");
	appendValue(seqs, "GGGGGGAGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");
	appendValue(seqs, "AAAAAATAGCTAGCTAGCTAGCTAGCTAGCTAGC");
	appendValue(seqs, "TAGTCATACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");

	StringSet<String<char> > ids;
	appendValue(ids, "Adenin1");
	appendValue(ids, "Guanin");
	appendValue(ids, "Adenin2");
	appendValue(ids, "Unidentifiziert");

	StringSet<String<Dna> > barcodes;
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "AAAAAA");

	DemultiplexStats demultiplexStats;
	resize(demultiplexStats.groups, length(barcodes)+1);

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

	Index<StringSet<String<Dna> >, IndexEsa<> > indexSet(barcodes);
	Finder<Index<StringSet<String<Dna> >, IndexEsa<> > > esaFinder(indexSet);
	std::vector<std::vector<int> > res = DoAll(seqs, barcodes, esaFinder, false, demultiplexStats);

	for(unsigned i = 0; i < length(exspected); ++i)
	{
		for(unsigned j = 0; j < length(exspected[i]); ++j)
		{
			SEQAN_ASSERT_EQ(exspected[i][j], res[i][j]);
			SEQAN_ASSERT_EQ(exspectedSeqs[exspected[i][j]], seqs[res[i][j]]);
		}
	}
}
// Checks the correctness of the DoAll function performing all demultiplexing operations for approximate inline barcode matching.
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
	
	DemultiplexStats demultiplexStats;
	resize(demultiplexStats.groups, length(barcodes)+1);

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

	std::vector<std::vector<int> > res = DoAll(seqs, barcodes, false, demultiplexStats);

	for(unsigned i = 0; i < length(exspected); ++i)
	{
		for(unsigned j = 0; j < length(exspected[i]); ++j)
		{
			SEQAN_ASSERT_EQ(exspected[i][j], res[i][j]);
			SEQAN_ASSERT_EQ(exspectedSeqs[exspected[i][j]], seqs[res[i][j]]);
		}
	}
}
// Checks the correctness of the DoAll function performing all demultiplexing operations for exact multiplex barcode matching.
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
	
	StringSet<String<Dna5Q> > barcodes;
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "AAAAAA");

	DemultiplexStats demultiplexStats;
	resize(demultiplexStats.groups, length(barcodes)+1);

	StringSet<String<Dna5Q> > exspectedSeqs;
	appendValue(exspectedSeqs, "GTACGATCGTACGTACGATGCTACGATGCATGCTACGATGCTACG");
	appendValue(exspectedSeqs, "AGTACGTACGTAGCTAGCTAGCATGCTAGCTAGCTAC");
	appendValue(exspectedSeqs, "TAGCTAGCTAGCTAGCTAGCTAGCTAGC");

	std::vector<std::vector<int> > exspected;
	resize(exspected, 4);
	appendValue(exspected[0], 1);	
	appendValue(exspected[3], 0);
	appendValue(exspected[3], 2);

	std::vector<std::vector<int> > res = DoAll(multiplex, barcodes, demultiplexStats);

	for(unsigned i = 0; i < length(exspected); ++i)
	{
		for(unsigned j = 0; j < length(exspected[i]); ++j)
		{
			SEQAN_ASSERT_EQ(exspected[i][j], res[i][j]);
			SEQAN_ASSERT_EQ(exspectedSeqs[exspected[i][j]], seqs[res[i][j]]);
		}
	}
}
// Checks the correctness of the DoAll function performing all demultiplexing operations for appromxiate multiplex barcode matching.
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

	StringSet<String<Dna5Q> > barcodes;
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "AAAAAA");

	DemultiplexStats demultiplexStats;
	resize(demultiplexStats.groups, length(barcodes)+1);

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

	std::vector<std::vector<int> > res = DoAll(multiplex, barcodes, demultiplexStats);

	for(unsigned i = 0; i < length(exspected); ++i)
	{
		for(unsigned j = 0; j < length(exspected[i]); ++j)
		{
			SEQAN_ASSERT_EQ(exspected[i][j], res[i][j]);
			SEQAN_ASSERT_EQ(exspectedSeqs[exspected[i][j]], seqs[res[i][j]]);
		}
	}
}
// Checks the correctness of the buildSets function which processes the result of the DoAll function so that the original sequences are not needed any longer. Also writes the sequences to a file.
SEQAN_DEFINE_TEST(Output_test)
{
	StringSet<String<Dna5Q> > seqs;
	appendValue(seqs, "AAAAAAGTACGATCGAAAAAAAAAAAAAAAAAAAAATGCTACGATGCTACG");
	appendValue(seqs, "GGGGGGAGTACGTACGGGGGGGGGGGGGAGCTAGCTAC");
	appendValue(seqs, "AAAAAATAGCTAGCTAGCTAAAAAAAAAAAAAAAGCTAGC");

	appendValue(seqs, "TAGTCATACGTACGTAGCNNNNNNNNNNNNNNNGCTAGCTAGCTAC");

	StringSet<String<char> > ids;
	appendValue(ids, "Adenin1");
	appendValue(ids, "Guanin");
	appendValue(ids, "Adenin2");
	appendValue(ids, "Unidentifiziert");

	StringSet<String<char> > barcodes;
	appendValue(barcodes, "GGGGGG");
	appendValue(barcodes, "CCCCCC");
	appendValue(barcodes, "AAAAAA");

	DemultiplexStats demultiplexStats;
	resize(demultiplexStats.groups, length(barcodes)+1);

	StringSet<String<char> > barcodeIds;
	appendValue(barcodeIds, "G-Barcode");
	appendValue(barcodeIds, "C-Barcode");
	appendValue(barcodeIds, "A-Barcode");

	StringSet<String<Dna5Q> > exspectedSeqs;
	appendValue(exspectedSeqs, "GTACGATCGAAAAAAAAAAAAAAAAAAAAATGCTACGATGCTACG");
	appendValue(exspectedSeqs, "AGTACGTACGGGGGGGGGGGGGAGCTAGCTAC");
	appendValue(exspectedSeqs, "TAGCTAGCTAGCTAAAAAAAAAAAAAAAGCTAGC");
	appendValue(exspectedSeqs, "TAGTCATACGTACGTAGCNNNNNNNNNNNNNNNGCTAGCTAGCTAC");

	std::vector<std::vector<int> > groups = DoAll(seqs, barcodes, false, demultiplexStats);
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
// Checks the correctness of the functions if they are applied to external data.
SEQAN_DEFINE_TEST(IO_test)
{
	StringSet<String<char> > ids; 
	StringSet<String<Dna5Q> > seqs;
	StringSet<String<char> > bcids;
	StringSet<String<Dna> > bcs;
	String<char> seqpath = "D:/seqan/test/seqs.fasta";
	String<char> bcpath = "D:/seqan/test/barcodes.fasta";

	SEQAN_ASSERT_EQ(0, loadSeqs(toCString(seqpath), ids, seqs));
	SEQAN_ASSERT_EQ(0, loadBarcodes(toCString(bcpath), bcids, bcs));

	DemultiplexStats demultiplexStats;
	resize(demultiplexStats.groups, length(bcs)+1);

	Index<StringSet<String<Dna> >, IndexEsa<> > indexSet(bcs);
	Finder<Index<StringSet<String<Dna> >, IndexEsa<> > > esaFinder(indexSet);
	std::vector<std::vector<int> > groups = DoAll(seqs, bcs, esaFinder, false, demultiplexStats);
	std::vector<StringSet<String<Dna5Q> > > gSeqs;
	std::vector<StringSet<String<char> > > gIds;

	buildSets(seqs, ids, groups, gSeqs, gIds);
	String<char> outpath =  "D:/seqan/out/iotest/";
	SEQAN_ASSERT_EQ(0, writeGroups(gSeqs, gIds, bcids, groups, outpath));
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
	SEQAN_CALL_TEST(IO_test);
}
SEQAN_END_TESTSUITE
