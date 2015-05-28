#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/file.h>
#include "adapterTrimming.h"

#define READFILE1 "/home/benni/Programmierung/Seqan/Daten/Sample_L3_Sample1/L3_Sample1_CGATGT_L003_R1_001.fastq.gz"
#define READFILE2 "/home/benni/Programmierung/Seqan/Daten/Sample_L3_Sample1/L3_Sample1_CGATGT_L003_R2_001.fastq.gz"

// Reads generated with SimSeq, so we know how long the inserts really are.
#define SIM_READFILE1 "/home/benni/Programmierung/Seqan/Daten/barcode.library.1.fastq"
#define SIM_READFILE2 "/home/benni/Programmierung/Seqan/Daten/barcode.library.2.fastq"

SEQAN_DEFINE_TEST(get_overlap_test)
{
	typedef seqan::String<seqan::Dna5> TSequence;
	typedef seqan::Align<TSequence, seqan::ArrayGaps > TAlign;
	typedef seqan::Row<TAlign>::Type TRow;

    TSequence seq1 = "AAAAAAAAA";
    TSequence seq2 = "TTTTTTT";

    TAlign align;
    resize(rows(align), 2);
    seqan::assignSource(row(align,0),seq1);
    seqan::assignSource(row(align,1),seq2);

    TRow &row1 = row(align, 0);
    TRow &row2 = row(align, 1);

    // Overlap of 4.
    // ---AAAAAAAAA
    // TTTTTTT-----
    seqan::insertGaps(row1, 0, 3);
    seqan::insertGaps(row2, length(row2), 5);

    SEQAN_ASSERT_EQ(getOverlap(align), 4u);

    // Overlap of 0.
    // -------AAAAAAAAA
    // TTTTTTT---------
    seqan::clearGaps(align);
    seqan::insertGaps(row1, 0, 7);
    seqan::insertGaps(row2, length(row2), 9);

    SEQAN_ASSERT_EQ(getOverlap(align), 0u);

    // Overlap of 7.
    // AAAAAAAAA
    // -TTTTTTT-
    seqan::clearGaps(align);
    seqan::insertGap(row2, 0);
    seqan::insertGap(row2, length(row2));

    SEQAN_ASSERT_EQ(getOverlap(align), 7u);
}

SEQAN_DEFINE_TEST(count_gap_test)
{
	seqan::Dna5String seq("AAGTCTATCTA");
	seqan::Gaps<seqan::Dna5String> row(seq);

	// no gaps yet.
	SEQAN_ASSERT_EQ(countTotalGaps(row), 0u);

	// insert a few gaps at the start.
	seqan::insertGaps(row, 0, 10);
	SEQAN_ASSERT_EQ(countTotalGaps(row), 10u);

	// create multiple continuous gaps

	seqan::insertGaps(row, 15, 10);
	seqan::insertGaps(row, length(row), 10);
	SEQAN_ASSERT_EQ(countTotalGaps(row), 30u);
}

SEQAN_DEFINE_TEST(insert_size_test)
{
	// The insert method is specified only for actually
	// overlapping sequences so we test those.
	typedef seqan::String<seqan::Dna5> TSequence;
	typedef seqan::Align<TSequence, seqan::ArrayGaps > TAlign;
	typedef seqan::Row<TAlign>::Type TRow;

    TSequence seq1 = "AAAAA";
    TSequence seq2 = "TTT";

    TAlign align;
    resize(rows(align), 2);
    seqan::assignSource(row(align,0),seq1);
    seqan::assignSource(row(align,1),seq2);

    TRow &row1 = row(align, 0);
    TRow &row2 = row(align, 1);

    // Insert size is 7
    // AAAAA--
    // ----TTT
    seqan::insertGaps(row1, seqan::length(row1), 2);
    seqan::insertGaps(row2, 0, 4);
    SEQAN_ASSERT_EQ(getInsertSize(align), 7u);

    // Insert size is 5
    // AAAAA
    // --TTT
    seqan::clearGaps(align);
    seqan::insertGaps(row2, 0, 2);
    SEQAN_ASSERT_EQ(getInsertSize(align), 5u);

    // Insert size is 4 (seq1 goes one base into the adapter).
	// AAAA|A
	// -TTT|-
	seqan::clearGaps(align);
	seqan::insertGap(row2, seqan::length(row2));
	seqan::insertGap(row2, 0);
	SEQAN_ASSERT_EQ(getInsertSize(align), 4u);
}

// Returns (num2 - num1 + 1) in a string "xxxx_num1_num2_xxx...".
int insertSize(char* s)
{
	int insert = 0;
	strtok(s, "_");
	insert -= atoi(strtok(NULL, "_"));
	insert += atoi(strtok(NULL, "_")) + 1;
	return insert;
}

SEQAN_DEFINE_TEST(simulated_insert_size_test)
{
	seqan::SequenceStream inStream1(SIM_READFILE1);
	seqan::SequenceStream inStream2(SIM_READFILE2);

	seqan::String<char> id1, id2;
	seqan::String<seqan::Dna5Q> seq1, seq2;

	typedef seqan::Align<seqan::String<seqan::Dna5Q> > TAlign;
	TAlign align;

	int count = 0, correct = 0;
	while (!(atEnd(inStream1) || atEnd(inStream2)))
	{
		if ((seqan::readRecord(id1, seq1, inStream1) +
			 seqan::readRecord(id2, seq2, inStream2)) == 0)
			{
				// Get actual insert size from simulated data.
				unsigned insert = insertSize(seqan::toCString(id1));
				// Ignore reads that don't overlap.
				if (insert > length(seq1) + length(seq2))
					continue;

				// Check if actual overlap equals estimated overlap.
				STRING_REVERSE_COMPLEMENT<seqan::Dna5Q>::Type mod(seq2);

				// Same alignment method stripPair uses.
				seqan::Pair<int, TAlign> ret = alignPair(seq1, mod, seqan::AlignConfig<true, true, true, true>());
				align = ret.i2;

				count++;
				if (insert == getInsertSize(align))
					correct++;
			}
	}

	// We should at least get 95% of all insert sizes right.
	// (The fact that we can't detect all comes from the overlap alignment.)
	double sensitvity = (double)correct/(double)count;
	// Create c-String out of double so we can print it in case of a problem.
	std::stringstream str;
	str << std::setprecision(2) << sensitvity;
	std::string tmp = str.str();
	SEQAN_ASSERT_MSG(sensitvity > 0.95, tmp.c_str());
}

SEQAN_DEFINE_TEST(match_test)
{
	Auto a = Auto();
	User u(7, 2);

	// Up to 5 overlap no error is allowed.
	SEQAN_ASSERT(isMatch(5,0,a));
	SEQAN_ASSERT_NOT(isMatch(5,1,a));
	// From 5 to 10 one error is allowed.
	SEQAN_ASSERT(isMatch(10,1,a));
	SEQAN_ASSERT_NOT(isMatch(10,2,a));
	// Otherwise 33% errors are allowed.
	SEQAN_ASSERT(isMatch(100,33,a));
	SEQAN_ASSERT_NOT(isMatch(100,34,a));

	// We need an overlap of length 7 and no more than 2 errors.
	SEQAN_ASSERT_NOT(isMatch(5,3,u));
	SEQAN_ASSERT_NOT(isMatch(7,3,u));
	SEQAN_ASSERT(isMatch(7,2,u));
}



SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
	SEQAN_CALL_TEST(get_overlap_test);
	SEQAN_CALL_TEST(count_gap_test);
	SEQAN_CALL_TEST(insert_size_test);
	SEQAN_CALL_TEST(match_test);
	//SEQAN_CALL_TEST(simulated_adapter_stripping_test);
	//SEQAN_CALL_TEST(simulated_paired_stripping_test);
	//SEQAN_CALL_TEST(simulated_paired_stripping_test2);
	SEQAN_CALL_TEST(simulated_insert_size_test);
}
SEQAN_END_TESTSUITE
