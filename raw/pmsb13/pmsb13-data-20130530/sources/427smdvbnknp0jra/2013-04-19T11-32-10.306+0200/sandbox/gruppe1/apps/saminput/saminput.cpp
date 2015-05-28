#include <iostream>
#include <seqan/bam_io.h>
#include <vector>

template <typename T>
int median(T& input)
{	
	std::nth_element( begin(input), begin(input) + length(input) / 2, end(input) ); 
	return *( begin(input) + length(input) / 2 ); 
}

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
		
	while (!atEnd(SamInput))
	{	
		readRecord(record, SamInput);
		std::cout << length(record.seq) << std::endl;
		appendValue(save,record.seq);
	}

	for (unsigned i=0;i<length(save);i++)
	{
		std::cout << i+1 << ". " << median(save) << "\n" << std::endl;
	}

 
    return 0;
}