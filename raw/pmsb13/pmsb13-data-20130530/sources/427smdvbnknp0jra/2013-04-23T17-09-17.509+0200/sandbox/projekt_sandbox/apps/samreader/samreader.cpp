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
 

   int l=0;
    for (unsigned i = 0; i < length(SamInput.header.sequenceInfos); ++i)
        l=l+(unsigned)SamInput.header.sequenceInfos[i].i1 << '\t' << SamInput.header.sequenceInfos[i].i2 << '\n';


    return 0;
}