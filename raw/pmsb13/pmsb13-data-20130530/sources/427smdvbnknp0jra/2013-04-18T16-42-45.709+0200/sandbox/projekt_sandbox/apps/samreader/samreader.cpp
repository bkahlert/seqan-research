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
        std::cout << bamInput.header.sequenceInfos[i].i1 << '\t' << bamInput.header.sequenceInfos[i].i2 << '\n';
		std::cout << "Test\n";
    // Open output stream, "-" means stdin on if reading, else stdout.
    seqan::BamStream bamOut("-", seqan::BamStream::WRITE);
    // Copy header.  The header is automatically written out before
    // the first record./*
    bamOut.header = bamInput.header;
    seqan::BamAlignmentRecord record;
	/*
	seqan::CharString value; 
	seqan::CharString tagValue;
	seqan::BamAlignmentRecord getTagValue(tagValue, "SN", record);
	*/
	
	seqan::StringSet<seqan::Dna5String> references;
	seqan::BamAlignmentRecord record;
// Read references and record.
	seqan::Align<seqan::Dna5String> align;
	if (record.rId != seqan::BamAlignmentRecord::INVALID_REFID)
    bamRecordToAlignment(align, references[record.rId], record);
/*
	unsigned lenSum = 0;
	while (!atEnd(bamInput))
    if (hasFlagUnmapped(record))
        lenSum += length(record.seq);
	*/
	
	//getTagType();
	//std::cout << getAlignmentLengthInRef(record);


	
 
    return 0;
}