#include <iostream>
#include <seqan/bam_io.h>

int main()
{
    // BamStream open input stream
    seqan::BamStream SamInput("/Informatik/Development/example.sam");
    if (!isGood(SamInput))
    {
		std::cerr << "-ERROR- Could not open" << std::endl;
		return 1;
    }
	// record represent a record from a SAM-file
	seqan::BamAlignmentRecord record;
	seqan::String<seqan::String<char>> save;
	
	unsigned count = 0;
	while (!atEnd(SamInput))
	{	
		readRecord(record, SamInput);
		std::cout << length(record.seq) << std::endl;
		appendValue(save,record.seq);
		count++;
	}
	std::cout << save << std::endl;

 
    return 0;
}