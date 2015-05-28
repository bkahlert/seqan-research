#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/file.h>
#include "adapterTrimming.h"

//#define READFILE1 "/home/benni/Programmierung/Seqan/Daten/Sample_L3_Sample1/L3_Sample1_CGATGT_L003_R1_001.fastq.gz"
//#define READFILE2 "/home/benni/Programmierung/Seqan/Daten/Sample_L3_Sample1/L3_Sample1_CGATGT_L003_R2_001.fastq.gz"
#define READFILE1 "/home/benni/Programmierung/Seqan/Daten/library.1.fastq"
#define READFILE2 "/home/benni/Programmierung/Seqan/Daten/library.2.fastq"

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
			int ret = stripAdapter(seq, adapter);
			trimmed_sum += ret;
			std::cout << "Davor: \n" << seq << "\n\n";
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
			alignPair(seq1, seq2);
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



SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    //SEQAN_CALL_TEST(adapter_alignment_test);
    //SEQAN_CALL_TEST(adapter_stripping_test);
	SEQAN_CALL_TEST(strip_pair_test);
}
SEQAN_END_TESTSUITE
;
