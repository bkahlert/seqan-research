#define SEQAN_PROFILE
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/file.h>

using namespace seqan;

int loadBatch(char* file, int records)
{
	int count = 0;
	SequenceStream inStream(file);
	if (!isGood(inStream))
	{
		std::cerr << "Could not open file.\n";
		return 1;
	}

	StringSet<CharString> ids;
	StringSet<String<Dna5Q> > seqs;

	while (!atEnd(inStream))
	{
		int status = readBatch(ids, seqs, inStream, records);
		count += length(seqs);
		if (status != 0)
		{
			std::cerr << "Error while reading file.\n";
			return 1;
		}
	}


	std::cout << "Read " << count << " sequences.\n";
	return 0;
}

int loadMultiSeq(char* file, int records)
{
	MultiSeqFile multiSeqFile;
	if (!open(multiSeqFile.concat, file, OPEN_RDONLY))
	{
		std::cerr << "Could not open file.\n";
		return 1;
	}

    AutoSeqFormat format;
    guessFormat(multiSeqFile.concat, format);
    split(multiSeqFile, format);

    unsigned seqCount = length(multiSeqFile);
    std::cout << seqCount << "\n";
    StringSet<String<Dna5Q> > seqs;
    StringSet<CharString> seqIDs;

    reserve(seqs, records, Exact());
    reserve(seqIDs, records, Exact());

    String<Dna5Q> seq;
    CharString qual;
    CharString id;

    int count = 0;
    for (unsigned r=0; r < seqCount; r += records)
    {
    	clear(seqs);
    	clear(seqIDs);
        reserve(seqs, records, Exact());
        reserve(seqIDs, records, Exact());

    	unsigned lim = _min(r+records, seqCount); // Prevent going out of bounds.
		for (unsigned i = r; i < lim; ++i)
		{
			++count;
			assignSeq(seq, multiSeqFile[i], format);    // read sequence
			assignQual(qual, multiSeqFile[i], format);  // read ascii quality values
			assignSeqId(id, multiSeqFile[i], format);   // read sequence id

			// convert ascii to values from 0..62
			// store dna and quality together in Dna5Q
			for (unsigned j = 0; j < length(qual) && j < length(seq); ++j)
				assignQualityValue(seq[j], (int)(ordValue(qual[j]) - 33));

			// we use reserve and append, as assign is not supported
			// by StringSet<..., Owner<ConcatDirect<> > >
			appendValue(seqs, seq, Generous());
			appendValue(seqIDs, id, Generous());
			//std::cout << seq << "\n";
		}
    }

    std::cout << "Read " << count << " sequences.\n";
    return 0;
}

int main(int argc, char** argv){
	if (argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " FILE.\n";
		return 1;
	}

	int records = 10000;

	SEQAN_PROTIMESTART(loadTime1);
	loadBatch(argv[1], records);
	std::cout << "Loading sequences took " << SEQAN_PROTIMEDIFF(loadTime1) << " seconds.\n";

	SEQAN_PROTIMESTART(loadTime2);
	loadMultiSeq(argv[1], records);
	std::cout << "Loading sequences took " << SEQAN_PROTIMEDIFF(loadTime2) << " seconds.\n";

	return 0;
}
