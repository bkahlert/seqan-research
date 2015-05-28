/*! Program part for managing input and output */
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

int loadSeqs(char const * path, StringSet<String<char> >& ids, StringSet<String<Dna5Q> >& seqs, StringSet<String<char> >& quals)
{
	SequenceStream seqStream(path, SequenceStream::READ);
	unsigned records = 1000;									//number of records to be read
	resize(ids, records);
	resize(seqs, records);
	resize(quals, records);

	if (!isGood(seqStream))
	{
		std::cerr << "Error while opening the sequence-file.\n";
		return 1;
	}
	
	if(readBatch(ids, seqs, quals, seqStream, records) != 0)		
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

int loadAdapters(char const * path, StringSet<String<char> >& adids, StringSet<String<Dna5> >& ads)
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

int loadMultiplexchar(char const * path, StringSet<String<Dna5Q> >& multiplex)
{
	SequenceStream seqStream(path, SequenceStream::READ);
	unsigned records = 1000;									//number of records to be read
	StringSet<String<char> > ids;
	StringSet<String<char> > quals;
	resize(multiplex, records);
	resize(ids, records);
	resize(quals, records);

	if (!isGood(seqStream))
	{
		std::cerr << "Error while opening the multiplex barcode-file.\n";
		return 1;
	}
	if(readBatch(ids, multiplex, quals, seqStream, records) != 0)		
	{
		std::cerr << "Error while reading the multiplex barcodes.\n";
		return 1;
	}
	return 0;
}

int writeSeqs(char const * path, const StringSet<String<char> >& ids, const StringSet<String<Dna5Q> >& seqs, const StringSet<String<char> >& quals)
{
	SequenceStream seqStream(path, SequenceStream::WRITE);
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not create/open the file.\n";
        return 1;
    }
	if(writeAll(seqStream, ids, seqs, quals) != 0)		
	{
		std::cerr << "Error while saving the sequences.\n";
		return 1;
	}
	return 0;
}