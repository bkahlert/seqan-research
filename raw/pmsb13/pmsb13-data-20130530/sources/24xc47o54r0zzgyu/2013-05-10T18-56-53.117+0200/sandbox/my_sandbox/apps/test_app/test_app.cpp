/*! Program part for parsing the input */
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

int loadSeqs(char const * path, StringSet<CharString>& ids, StringSet<Dna5String>& seqs, StringSet<CharString>& quals)
{
	SequenceStream seqStream(path, SequenceStream::READ);
	unsigned reads = 1000;									//number of records to be read
	resize(ids, reads);
	resize(seqs, reads);
	resize(quals, reads);

	if (!isGood(seqStream))
	{
		std::cerr << "Error while opening the sequence-file.\n";
		return 1;
	}
	
	if(readBatch(ids, seqs, quals, seqStream, reads) != 0)		
	{
		std::cerr << "Error while reading the sequences.\n";
		return 1;
	}
	return 0;
}

int loadBarcodes(char const * path, StringSet<CharString>& bcids, StringSet<DnaString>& bcs)
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

int loadAdapters(char const * path, StringSet<CharString>& adids, StringSet<DnaString>& ads)
{
	SequenceStream adStream(path, SequenceStream::READ);

	if(!isGood (adStream))
	{
		std::cerr << "Error while opening adapter-file.\n";
		return 1;
	}

	if(readAll(adids, ads, adStream) != 0)
	{
		std::cerr << "Error while reading the adapters.\n";
		return 1;
	}
	return 0;
}

int main(int argc, char const * argv[])
{
	StringSet<CharString> bcids;
	StringSet<DnaString> bcs;
	
	loadBarcodes(argv[1], bcids, bcs);
	for (unsigned i = 0; i < 5; ++i)
	{
		std::cout << bcids[i] << "\t" << bcs[i] <<std::endl;
	}
}