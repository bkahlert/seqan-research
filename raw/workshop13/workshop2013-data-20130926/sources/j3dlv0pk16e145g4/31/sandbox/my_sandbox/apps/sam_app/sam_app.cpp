#include <iostream>
#include <seqan/bam_io.h>

int main()
{
    // Open input stream, BamStream can read SAM and BAM files.
    seqan::BamStream bamStreamIn("/home/felix/test.sam");
    if (!isGood(bamStreamIn))
    {
        std::cerr << "ERROR: Could not open test.sam.";
        return 1;
    }
    // Open output stream, "-" means stdin on if reading, else stdout.
    seqan::BamStream bamStreamOut("-", seqan::BamStream::WRITE);
    // Copy header.  The header is automatically written out before
    // the first record.
    bamStreamOut.header = bamStreamIn.header;
    seqan::BamAlignmentRecord record;
    while (!atEnd(bamStreamIn))
    {
        readRecord(record, bamStreamIn);
        writeRecord(bamStreamOut, record);
    }
    return 0;
}