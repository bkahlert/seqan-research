#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;
using namespace std;

int main(int argc, char** argv )
{
    if (argc>2)
	return 1;
    
    
    seqan::SequenceStream seqStream(argv[1], seqan::SequenceStream::WRITE);
    
    if (!isGood(seqStream))
    {
	std::cerr << "ERROR: Could not open the file.\n";
	return 1;
    }
    
    
    seqan::String<CharString> ids;
    seqan::String<Dna5String> seqs;
    seqan::String<CharString> quals;
    
    typedef Iterator<String<CharString> >::Type TCharStringIterator;
    typedef Iterator<String<Dna5String> >::Type TDna5StringIterator;
    
    TCharStringIterator it1=begin(ids);
    TDna5StringIterator it2=begin(seqs);
    for (;((!atEnd(it1))&&(!atEnd(it2)));goNext(it1),goNext(it2))
	if (writeRecord(seqStream, *it1, *it2) != 0)
	{
	    std::cerr << "ERROR: Could not write to " << argv[1] << "!\n";
	    return 1;
	}
    
    return 0;
}