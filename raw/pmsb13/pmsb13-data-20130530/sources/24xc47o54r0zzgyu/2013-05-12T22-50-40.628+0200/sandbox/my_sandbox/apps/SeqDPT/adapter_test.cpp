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
#define SIM_READFILE1 "/home/benni/Programmierung/Seqan/Daten/library.1.fastq"
#define SIM_READFILE2 "/home/benni/Programmierung/Seqan/Daten/library.2.fastq"

// No actual testing yet, only print alignments for inspection.
SEQAN_DEFINE_TEST(adapter_alignment_test)
{
	seqan::String<seqan::Dna5> adapter("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG");

	// No error checking, we assume the file exists (it's a constant test file).
	seqan::SequenceStream inStream(READFILE1);
	seqan::String<char> id;
	seqan::String<seqan::Dna5> seq;
	seqan::CharString qual;

	int i=0;
	while (!atEnd(inStream) && i++ < 200)
	{
		if (seqan::readRecord(id, seq, qual, inStream) == 0)
		{
			int score = alignAdapter(seq, adapter);
			if (score > 5)
				std::cout << score << "\n";
		}
	}
}

SEQAN_DEFINE_TEST(adapter_stripping_test)
{
	seqan::String<seqan::Dna5> adapter("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG");
	//seqan::String<seqan::Dna5> adapter("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT");

	// No error checking, we assume the file exists (it's a constant test file).
	seqan::SequenceStream inStream(READFILE1);

	if (!seqan::isGood(inStream))
	{
		std::cout << "Fehler beim Öffnen von Datei" << std::endl;
	}

	seqan::String<char> id;
	seqan::String<seqan::Dna5> seq;
	seqan::CharString qual;
	int i=0;
	int count = 0, sum = 0, trimmed_sum = 0;
	while (!atEnd(inStream) && i++ < 200)
	{
		if (seqan::readRecord(id, seq, qual, inStream) == 0)
		{
			sum += length(seq);
			Auto a = Auto();
			std::cout << "Davor: \n" << seq << "\n\n";
			int ret = stripAdapter(seq, adapter, a);
			trimmed_sum += ret;
			if (ret > 0)
			{
				std::cout << "Danach(" << ret << "):\n"<< seq << "\n";
				++count;
			}
		}

		if (i % 10000 == 0)
			std::cout << i << "| Getrimmte Adapter: " << count << ", Anteil getrimmt: " << (double)trimmed_sum/(double)sum << "\n";
	}
}

SEQAN_DEFINE_TEST(align_pair_test)
{
	seqan::SequenceStream inStream1(READFILE1);
	seqan::SequenceStream inStream2(READFILE2);

	seqan::String<char> id1, id2;
	seqan::String<seqan::Dna5> seq1, seq2;
	seqan::CharString qual1, qual2;

	int i=0;
	while (!(atEnd(inStream1) || atEnd(inStream2)) && i++ < 200)
	{
		if ((seqan::readRecord(id1, seq1, qual1, inStream1) +
			seqan::readRecord(id2, seq2, qual2, inStream2)) == 0)
		{
			alignPair(seq1, seq2, seqan::AlignConfig<true, true, true, true>());
		}
	}
}

SEQAN_DEFINE_TEST(strip_pair_test)
{
	seqan::SequenceStream inStream1(READFILE1);
	seqan::SequenceStream inStream2(READFILE2);

	seqan::String<char> id1, id2;
	seqan::String<seqan::Dna5> seq1, seq2;
	seqan::CharString qual1, qual2;

	int i=0;
	while (!(atEnd(inStream1) || atEnd(inStream2)) && i++ < 200)
	{
		if ((seqan::readRecord(id1, seq1, qual1, inStream1) +
			seqan::readRecord(id2, seq2, qual2, inStream2)) == 0)
		{
			std::cout << seq1 << "\n" << seq2 << "\n SCHNITT \n";
			stripPair(seq1, seq2);
			std::cout << seq1 << "\n" << seq2 << "\n\n" << i << "\n";
		}
	}
}

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

    SEQAN_ASSERT_EQ(getOverlap(align), 4);

    // Overlap of 0.
    // -------AAAAAAAAA
    // TTTTTTT---------
    seqan::clearGaps(align);
    seqan::insertGaps(row1, 0, 7);
    seqan::insertGaps(row2, length(row2), 9);

    SEQAN_ASSERT_EQ(getOverlap(align), 0);

    // Overlap of 7.
    // AAAAAAAAA
    // -TTTTTTT-
    seqan::clearGaps(align);
    seqan::insertGap(row2, 0);
    seqan::insertGap(row2, length(row2));

    SEQAN_ASSERT_EQ(getOverlap(align), 7);
}

SEQAN_DEFINE_TEST(count_gap_test)
{
	seqan::Dna5String seq("AAGTCTATCTA");
	seqan::Gaps<seqan::Dna5String> row(seq);

	// no gaps yet.
	SEQAN_ASSERT_EQ(countTotalGaps(row), 0);

	// insert a few gaps at the start.
	seqan::insertGaps(row, 0, 10);
	SEQAN_ASSERT_EQ(countTotalGaps(row), 10);

	// create multiple continuous gaps

	seqan::insertGaps(row, 15, 10);
	seqan::insertGaps(row, length(row), 10);
	SEQAN_ASSERT_EQ(countTotalGaps(row), 30);
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
    seqan::insertGaps(row1, length(row1), 2);
    seqan::insertGaps(row2, 0, 4);
    SEQAN_ASSERT_EQ(getInsertSize(align), 7);

    // Insert size is 5
    // AAAAA
    // --TTT
    seqan::clearGaps(align);
    seqan::insertGaps(row2, 0, 2);
    SEQAN_ASSERT_EQ(getInsertSize(align), 5);

    // Insert size is 4 (seq1 goes one base into the adapter).
	// AAAA|A
	// -TTT|-
	seqan::clearGaps(align);
	seqan::insertGap(row2, length(row2));
	seqan::insertGap(row2, 0);
	SEQAN_ASSERT_EQ(getInsertSize(align), 4);

}

SEQAN_DEFINE_TEST(simulated_insert_size_test)
{
	seqan::SequenceStream inStream1(SIM_READFILE1);
	seqan::SequenceStream inStream2(SIM_READFILE2);

	seqan::String<char> id1, id2;
	seqan::String<seqan::Dna5Q> seq1, seq2;

	while (!(atEnd(inStream1) || atEnd(inStream2)))
	{
		if ((seqan::readRecord(id1, seq1, inStream1) +
			seqan::readRecord(id2, seq2, inStream2)) == 0)
			{
				// Get actual insert size from simulated data.

				// Ignore reads that don't overlap.

				// Check if actual overlap equals estimated overlap.
			}
	}
}


SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    //SEQAN_CALL_TEST(adapter_alignment_test);
    //SEQAN_CALL_TEST(adapter_stripping_test);
	//SEQAN_CALL_TEST(strip_pair_test);
	SEQAN_CALL_TEST(get_overlap_test);
	SEQAN_CALL_TEST(count_gap_test);
	SEQAN_CALL_TEST(insert_size_test);
}
SEQAN_END_TESTSUITE
;
