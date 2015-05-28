#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/file.h>
#include "readTrimming.h"


#define TESTSAMPLE "/home/benni/Programmierung/Seqan/Daten/testsample.fq"

SEQAN_DEFINE_TEST(sliding_window_test_qual)
{
	// No error checking, we assume the file exists (it's a constant test file).
	seqan::SequenceStream inStream(TESTSAMPLE);
	seqan::String<char> id;
	seqan::String<seqan::Dna5> seq;
	seqan::CharString qual;

	// The number of bases Trimmomatic cuts from the sample file with the same settings.
	unsigned expected[] = {12, 21, 7, 75, 100, 100, 4, 6, 3, 3, 3, 100, 3, 2, 57, 2, 3, 31, 97, 52,
					   2, 4, 14, 2, 2, 3, 2, 3, 4, 2, 70, 3, 2, 3, 3, 3, 3, 54, 2, 8, 3, 3, 98, 3};

	unsigned i=0;
	Mean m = Mean(5);
	while (!atEnd(inStream))
	{
		if (seqan::readRecord(id, seq, qual, inStream) == 0)
		{
			unsigned len = length(seq);
			unsigned res = trimRead(seq, qual, 20, m);
			SEQAN_ASSERT_EQ(expected[i++], res);     // Check if returned number of trimmed bases is correct.
			SEQAN_ASSERT_EQ(res, len - length(seq)); // Check if we really trimmed so many bases.
		}
	}
}

SEQAN_DEFINE_TEST(drop_reads_test)
{
	typedef seqan::String<seqan::Dna5Q> TString;
	typedef seqan::CharString TID;
	seqan::StringSet<TString> test;
	seqan::StringSet<TID> testID;

	seqan::appendValue(testID, TString("1"));
	seqan::appendValue(test, TString("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
	seqan::appendValue(testID, TString("2"));
	seqan::appendValue(test, TString("AAAAAAAAAAAAAAAAAAAA"));
	seqan::appendValue(testID, TString("3"));
	seqan::appendValue(test, TString("AAAAA"));
	seqan::appendValue(testID, TString("4"));
	seqan::appendValue(test, TString("AAA"));
	seqan::appendValue(testID, TString("5"));
	seqan::appendValue(test, TString("A"));

	QualityTrimmingStats tmp;

	dropReads(testID, test, 2, tmp);
	SEQAN_ASSERT_EQ(length(test), 4u);

	dropReads(testID, test, 4, tmp);
	SEQAN_ASSERT_EQ(length(test), 3u);

	dropReads(testID, test, 6, tmp);
	SEQAN_ASSERT_EQ(length(test), 2u);

	dropReads(testID, test, 24, tmp);
	SEQAN_ASSERT_EQ(length(test), 1u);

	dropReads(testID, test, 50, tmp);
	SEQAN_ASSERT_EQ(length(test), 0u);
}

SEQAN_DEFINE_TEST(remove_value_test)
{
	typedef seqan::String<seqan::Dna5Q> TString;
	typedef seqan::CharString TID;
	seqan::StringSet<TString> test;
	seqan::StringSet<TID> testID;

	seqan::appendValue(testID, TString("1"));
	seqan::appendValue(test, TString("ACGT"));
	seqan::appendValue(testID, TString("2"));
	seqan::appendValue(test, TString("TGCA"));

	// Remove the first string and its ID.
	removeValues(test, testID, 0);
	// Make sure we removed something.
	SEQAN_ASSERT_EQ(length(test), 1u);
	// Make sure we removed the correct values.
	SEQAN_ASSERT_EQ(seqan::value(test,0), TString("TGCA"));
}


SEQAN_DEFINE_TEST(sliding_window_test)
{
	// No error checking, we assume the file exists (it's a constant test file).
	seqan::SequenceStream inStream(TESTSAMPLE);
	seqan::String<char> id;
	seqan::String<seqan::Dna5Q> seq;

	// The number of bases Trimmomatic cuts from the sample file with the same settings.
	unsigned expected[] = {12, 21, 7, 75, 100, 100, 4, 6, 3, 3, 3, 100, 3, 2, 57, 2, 3, 31, 97, 52,
					   2, 4, 14, 2, 2, 3, 2, 3, 4, 2, 70, 3, 2, 3, 3, 3, 3, 54, 2, 8, 3, 3, 98, 3};

	unsigned i=0;
	Mean m = Mean(5);
	while (!atEnd(inStream))
	{
		if (seqan::readRecord(id, seq, inStream) == 0)
		{
			unsigned len = length(seq);
			unsigned res = trimRead(seq, 20, m);
			SEQAN_ASSERT_EQ(expected[i++], res);	 // Check if returned number of trimmed bases is correct.
			SEQAN_ASSERT_EQ(res, len - length(seq)); // Check if we really trimmed so many bases.
		}
	}
}

SEQAN_DEFINE_TEST(cut_tail_test)
{
	unsigned expected[] = {9, 9, 3, 26, 69, 3, 2, 0, 0, 0, 0, 23, 1, 0, 22, 0, 0, 12, 35, 6,
					   0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 52, 0};

	seqan::SequenceStream inStream(TESTSAMPLE);
	seqan::String<char> id;
	seqan::String<seqan::Dna5> seq;
	seqan::CharString qual;

	unsigned i=0;
	Tail t;
	while (!atEnd(inStream))
	{
		if (seqan::readRecord(id, seq, qual, inStream) == 0)
		{
			unsigned len = length(seq);
			unsigned res = trimRead(seq, qual, 20, t);
			SEQAN_ASSERT_EQ(expected[i++], res);     // Check if returned number of trimmed bases is correct.
			SEQAN_ASSERT_EQ(res, len - length(seq)); // Check if we really trimmed so many bases.
		}
	}
}

SEQAN_DEFINE_TEST(cut_bwa_test)
{
	unsigned expected[] = {8, 8, 2, 25, 68, 2, 1, 0, 0, 0, 0, 22, 0, 0, 21, 0, 0, 11, 34, 5,
					   0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 51, 0};

	seqan::SequenceStream inStream(TESTSAMPLE);
	seqan::String<char> id;
	seqan::String<seqan::Dna5> seq;
	seqan::CharString qual;

	unsigned i=0;
	BWA b;
	while (!atEnd(inStream))
	{
		if (seqan::readRecord(id, seq, qual, inStream) == 0)
		{
			unsigned len = length(seq);
			unsigned res = trimRead(seq, qual, 20, b);
			SEQAN_ASSERT_EQ(expected[i++], res);     // Check if returned number of trimmed bases is correct.
			SEQAN_ASSERT_EQ(res, len - length(seq)); // Check if we really trimmed so many bases.
		}
	}
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(sliding_window_test_qual);
    SEQAN_CALL_TEST(sliding_window_test);
    SEQAN_CALL_TEST(drop_reads_test);
    SEQAN_CALL_TEST(remove_value_test);
    SEQAN_CALL_TEST(cut_tail_test);
    SEQAN_CALL_TEST(cut_bwa_test);
}
SEQAN_END_TESTSUITE
