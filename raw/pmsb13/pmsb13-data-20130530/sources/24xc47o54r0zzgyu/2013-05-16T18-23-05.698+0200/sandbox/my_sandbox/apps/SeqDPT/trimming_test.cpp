#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/file.h>
#include "readTrimming.h"


SEQAN_DEFINE_TEST(sliding_window_test_qual)
{
	// No error checking, we assume the file exists (it's a constant test file).
	seqan::SequenceStream inStream("/home/benni/Programmierung/Seqan/Daten/testsample.fq");
	seqan::String<char> id;
	seqan::String<seqan::Dna5> seq;
	seqan::CharString qual;

	// The number of bases Trimmomatic cuts from the sample file with the same settings.
	unsigned expected[] = {12, 21, 7, 75, 100, 100, 4, 6, 3, 3, 3, 100, 3, 2, 57, 2, 3, 31, 97, 52,
					   2, 4, 14, 2, 2, 3, 2, 3, 4, 2, 70, 3, 2, 3, 3, 3, 3, 54, 2, 8, 3, 3, 98, 3};

	int i=0;
	Mean m = Mean(5);
	while (!atEnd(inStream))
	{
		if (seqan::readRecord(id, seq, qual, inStream) == 0)
		{
			unsigned res = trimRead(seq, qual, 20, m);
			SEQAN_ASSERT_EQ(expected[i++], res);
		}
	}
}

SEQAN_DEFINE_TEST(sliding_window_test)
{
	// No error checking, we assume the file exists (it's a constant test file).
	seqan::SequenceStream inStream("/home/benni/Programmierung/Seqan/Daten/testsample.fq");
	seqan::String<char> id;
	seqan::String<seqan::Dna5Q> seq;

	// The number of bases Trimmomatic cuts from the sample file with the same settings.
	int expected[] = {12, 21, 7, 75, 100, 100, 4, 6, 3, 3, 3, 100, 3, 2, 57, 2, 3, 31, 97, 52,
					   2, 4, 14, 2, 2, 3, 2, 3, 4, 2, 70, 3, 2, 3, 3, 3, 3, 54, 2, 8, 3, 3, 98, 3};

	int i=0;
	Mean m = Mean(5);
	while (!atEnd(inStream))
	{
		if (seqan::readRecord(id, seq, inStream) == 0)
		{
			int res = trimRead(seq, 20, m);
			SEQAN_ASSERT_EQ(expected[i++], res);
		}
	}
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(sliding_window_test_qual);
    SEQAN_CALL_TEST(sliding_window_test);
}
SEQAN_END_TESTSUITE
