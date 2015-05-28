#include <iostream>
#include <seqan/bam_io.h>

template <typename T>
T median(T& input)
{	
	std::nth_element( begin(input), begin(input) + length(input) / 2, end(input) ); 
	return *( begin(input) + length(input) / 2 ); 
}

template <typename T>
int average(T input)
{	
	int tmpi=0;
	for (unsigned i=0;i<length(input);i++)
	{
		tmpi = tmpi + length(input[i]);
	}
	return tmpi/length(input);
}

int saminput(seqan::String<seqan::String<char>> &save, char *file)
{
    // BamStream open input stream
    seqan::BamStream SamInput(file);
    if (!isGood(SamInput))
    {
		std::cerr << "-ERROR- Could not open" << std::endl;
		return 1;
    }
	// record represent a record from a SAM-file
	seqan::BamAlignmentRecord record;

	while (!atEnd(SamInput))
	{	
		readRecord(record, SamInput);
		std::cout << length(record.seq) << std::endl;
		appendValue(save,record.seq);		
	}

	for (unsigned i=0;i<length(save);i++)
	{
		std::cout << i+1 << ". " << save[i] << "\n" << std::endl;
	}

 
    return 0;
}
template <typename T>
T devide(T input)
{
	int tmpi = average(input);
	T less;
	T great;
	for (unsigned i=0;i<length(input);i++)
	{
		if(length(input[i])<tmpi){append(less,input[i]);}
		else{append(great,input[i]);}
	}
	append(less,great);
	return less;
}

int main(int argc, char *argv[])
{
	argv[1] = "/Informatik/Development/example.sam";
	seqan::String<seqan::String<char>> save;
	saminput(save,argv[1]);

	save=devide(save);

	
	for (unsigned i=0;i<length(save);i++)
	{
		std::cout << save[i] << "\n" << std::endl;
	}


	return 0;
}