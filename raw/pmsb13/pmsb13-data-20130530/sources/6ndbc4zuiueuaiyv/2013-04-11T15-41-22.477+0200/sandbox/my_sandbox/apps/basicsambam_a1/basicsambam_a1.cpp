#include <iostream>
#include <seqan/bam_io.h>

using namespace seqan;

int main()
{
    // Open input stream, BamStream can read SAM and BAM files.
    seqan::BamStream bamStreamIn("/home/stefan/Dokumente/subversion/seqan-trunk/sandbox/my_sandbox/apps/basicsambam_a1/example.sam");
    // Open output stream, "-" means stdin on if reading, else stdout.
    seqan::BamStream bamStreamOut("-", seqan::BamStream::WRITE);
    // Copy header.  The header is automatically written out before
    // the first record.
    bamStreamOut.header = bamStreamIn.header;

    seqan::BamAlignmentRecord record;
    while (!atEnd(bamStreamIn) && (readRecord(record, bamStreamIn))==0 && (writeRecord(bamStreamOut, record))==0)
    
    return 0;
}

