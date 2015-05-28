#include <iostream>
#include <seqan/bam_io.h>


int main()
{
    // Open input stream, BamStream can read SAM and BAM files.
    seqan::BamStream bamInput("/Informatik/Development/example.sam");
    if (!isGood(bamInput))
    {
		std::cerr << "-ERROR- Could not open" << std::endl;
		return 1;
    }
	//prints the sequences
	for (unsigned i = 0; i < length(bamInput.header.sequenceInfos); ++i)
        std::cout << bamInput.header.sequenceInfos[i].i1 << '\t'
                  << bamInput.header.sequenceInfos[i].i2 << '\n';
    // Open output stream, "-" means stdin on if reading, else stdout.
    seqan::BamStream bamOut("-", seqan::BamStream::WRITE);
    // Copy header.  The header is automatically written out before
    // the first record.
    bamOut.header = bamInput.header;
    seqan::BamAlignmentRecord record;

    while (!atEnd(bamInput))
    {
        if (readRecord(record, bamInput) != 0)
        {
            std::cerr << "-ERROR- Could not read record" << std::endl;
            return 1;
        }
        if (writeRecord(bamOut, record) != 0)
        {
            std::cerr << "-ERROR- Could not write record!" << std::endl;
            return 1;
        }
    }

    return 0;
}