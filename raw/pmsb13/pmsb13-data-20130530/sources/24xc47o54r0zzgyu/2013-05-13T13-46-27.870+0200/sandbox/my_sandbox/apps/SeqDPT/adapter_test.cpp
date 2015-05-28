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
	seqan::String<seqan::Dna5Q> seq1, seq2;

	int i=0;
	while (!(atEnd(inStream1) || atEnd(inStream2)) && i++ < 200)
	{
		if ((seqan::readRecord(id1, seq1, inStream1) +
			seqan::readRecord(id2, seq2, inStream2)) == 0)
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
				int insert = insertSize(seqan::toCString(id1));
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

SEQAN_DEFINE_TEST(simulated_adapter_stripping_test)
{
	seqan::SequenceStream inStream(SIM_READFILE1);
	seqan::String<seqan::Dna5Q> adapter("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG");
	seqan::String<char> id;
	seqan::String<seqan::Dna5Q> seq;

	Auto a = Auto();

	int true_pos = 0, true_neg = 0,
        false_pos = 0, false_neg = 0;
    int count = 0;

	while (!atEnd(inStream))
	{
		if (readRecord(id, seq, inStream) == 0)
		{
			++count;
			int insert = insertSize(seqan::toCString(id));
			// We know how much adapter sequence we have read.
			int overlap = length(seq) > insert ? length(seq) - insert : 0;
			int trimmed = stripAdapter(seq, adapter, a);

			if (trimmed > overlap)  // Trimmed more than necessary.
				++false_pos;
			else if (trimmed == overlap && overlap > 0)  // Trimmed the adapter exactly.
				++true_pos;
			else if (trimmed < overlap)  // Trimmed less than necessary (or nothing at all).
				++false_neg;
			else if(trimmed == 0 && overlap == 0) // Correctly trimmed nothing.
				++true_neg;
		}
	}

	double sensitivity = (double) true_pos / ((double) true_pos + (double) false_neg);
	double specifivity = (double) true_neg / ((double) true_neg + (double) false_pos);

	double exact = (double) (true_pos + true_neg) / (double) count;
	double wrong = (double) (false_pos + false_neg) / (double) count;

	std::cout << "Sens: " << sensitivity << " Spec: " << specifivity << "\n";
	std::cout << "Exact: " << exact << " Wrong: " << wrong << "\n";
}

SEQAN_DEFINE_TEST(simulated_paired_stripping_test)
{
	seqan::SequenceStream inStream1(SIM_READFILE1);
	seqan::SequenceStream inStream2(SIM_READFILE2);
	seqan::String<char> id1, id2;
	seqan::String<seqan::Dna5Q> seq1, seq2;

	int true_pos = 0, true_neg = 0,
        false_pos = 0, false_neg = 0;
    int count = 0;

	while (!(atEnd(inStream1) || atEnd(inStream2)))
	{
		if (readRecord(id1, seq1, inStream1) + readRecord(id2, seq2, inStream2) == 0)
		{
			++count;
			int insert = insertSize(seqan::toCString(id1));
			int len1 = length(seq1);
			int len2 = length(seq2);
			// We know how much adapter sequence we have read.
			int overlap1 = length(seq1) > insert ? length(seq1) - insert : 0;
			int overlap2 = length(seq2) > insert ? length(seq2) - insert : 0;
			int est_insert = stripPair(seq1, seq2);
			int trimmed = len1 - length(seq1);

			if (trimmed > overlap1)  // Trimmed more than necessary.
				++false_pos;
			else if (trimmed == overlap1 && overlap1 > 0)  // Trimmed the adapter exactly.
				++true_pos;
			else if (trimmed < overlap1)  // Trimmed less than necessary (or nothing at all).
				++false_neg;
			else if(trimmed == 0 && overlap1 == 0) // Correctly trimmed nothing.
				++true_neg;
		}
	}

	double sensitivity = (double) true_pos / ((double) true_pos + (double) false_neg);
	double specifivity = (double) true_neg / ((double) true_neg + (double) false_pos);

	double exact = (double) (true_pos + true_neg) / (double) count;
	double wrong = (double) (false_pos + false_neg) / (double) count;

	std::cout << "Sens: " << sensitivity << " Spec: " << specifivity << "\n";
	std::cout << "Exact: " << exact << " Wrong: " << wrong << "\n";

	SEQAN_ASSERT(sensitivity > 0.95);
	SEQAN_ASSERT(specifivity > 0.95);
}

SEQAN_DEFINE_TEST(simulated_paired_stripping_test2)
{
	seqan::String<seqan::Dna5Q> adapter1("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG");
	seqan::String<seqan::Dna5Q> adapter2("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT");

	seqan::SequenceStream inStream1(SIM_READFILE1);
	seqan::SequenceStream inStream2(SIM_READFILE2);
	seqan::String<char> id1, id2;
	seqan::String<seqan::Dna5Q> seq1, seq2;

	int true_pos = 0, true_neg = 0,
        false_pos = 0, false_neg = 0;
    int count = 0;

    int i = 0;
	while (!(atEnd(inStream1) || atEnd(inStream2)))// && i++ < 100)
	{
		if (readRecord(id1, seq1, inStream1) + readRecord(id2, seq2, inStream2) == 0)
		{
			++count;
			int insert = insertSize(seqan::toCString(id1));
			int len1 = length(seq1);
			int len2 = length(seq2);
			// We know how much adapter sequence we have read.
			int overlap1 = length(seq1) > insert ? length(seq1) - insert : 0;
			int overlap2 = length(seq2) > insert ? length(seq2) - insert : 0;
			int est_insert = stripPair(seq1, adapter1, seq2, adapter2);

			//std::cout << insert << " " << est_insert << " " << insert - len1 - len2 << "\n";

			int trimmed = len1 - length(seq1);

			if (trimmed > overlap1)  // Trimmed more than necessary.
				++false_pos;
			else if (trimmed == overlap1 && overlap1 > 0)  // Trimmed the adapter exactly.
				++true_pos;
			else if (trimmed < overlap1)  // Trimmed less than necessary (or nothing at all).
				++false_neg;
			else if(trimmed == 0 && overlap1 == 0) // Correctly trimmed nothing.
				++true_neg;
		}
	}

	double sensitivity = (double) true_pos / ((double) true_pos + (double) false_neg);
	double specifivity = (double) true_neg / ((double) true_neg + (double) false_pos);

	double exact = (double) (true_pos + true_neg) / (double) count;
	double wrong = (double) (false_pos + false_neg) / (double) count;

	std::cout << "Sens: " << sensitivity << " Spec: " << specifivity << "\n";
	std::cout << "Exact: " << exact << " Wrong: " << wrong << "\n";

	std::cout << true_pos << " " << true_neg << " "<< false_pos << " "<< false_neg << "\n";

	//SEQAN_ASSERT(sensitivity > 0.95);
	//SEQAN_ASSERT(specifivity > 0.95);
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
	// Otherwise 15% errors are allowed.
	SEQAN_ASSERT(isMatch(100,15,a));
	SEQAN_ASSERT_NOT(isMatch(100,16,a));

	// We need an overlap of length 7 and no more than 2 errors.
	SEQAN_ASSERT_NOT(isMatch(5,3,u));
	SEQAN_ASSERT_NOT(isMatch(7,3,u));
	SEQAN_ASSERT(isMatch(7,2,u));
}



SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    //SEQAN_CALL_TEST(adapter_alignment_test);
    //SEQAN_CALL_TEST(adapter_stripping_test);
	//SEQAN_CALL_TEST(strip_pair_test);
	SEQAN_CALL_TEST(get_overlap_test);
	SEQAN_CALL_TEST(count_gap_test);
	SEQAN_CALL_TEST(insert_size_test);
	SEQAN_CALL_TEST(match_test);
	//SEQAN_CALL_TEST(simulated_adapter_stripping_test);
	//SEQAN_CALL_TEST(simulated_paired_stripping_test);
	SEQAN_CALL_TEST(simulated_paired_stripping_test2);
	//SEQAN_CALL_TEST(simulated_insert_size_test);
}
SEQAN_END_TESTSUITE
;
