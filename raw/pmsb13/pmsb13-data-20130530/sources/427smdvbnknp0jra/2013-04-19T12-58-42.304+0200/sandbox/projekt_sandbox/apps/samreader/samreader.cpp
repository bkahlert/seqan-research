#include <iostream>
#include <seqan/bam_io.h>

int main()
{
    // Open input stream, BamStream can read SAM and BAM files.
    seqan::BamStream SamInput("/Informatik/Development/example.sam");
    if (!isGood(SamInput))
    {
		std::cerr << "-ERROR- Could not open" << std::endl;
		return 1;
    }
	seqan::BamAlignmentRecord record;
	while (!atEnd(SamInput))
	{
    readRecord(record, SamInput);
	std::cout << length(record.seq) << std::endl;
	}
 
	{
    seqan::BamStream bamInStream("example.bam");

    for (unsigned i = 0; i < length(bamInStream.header.sequenceInfos); ++i)
        std::cout << bamInStream.header.sequenceInfos[i].i1 << '\t'
                  << bamInStream.header.sequenceInfos[i].i2 << '\n';

    return 0;
}
    return 0;
}