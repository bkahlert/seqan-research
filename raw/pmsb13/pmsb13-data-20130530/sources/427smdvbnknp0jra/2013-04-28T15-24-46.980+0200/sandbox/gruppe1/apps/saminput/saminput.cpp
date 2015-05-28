#include <iostream>
#include <seqan/bam_io.h>
#include <math.h>

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
		tmpi += length(input[i]);
	}
	return tmpi/length(input);
}

int analyse(seqan::String<seqan::String<char>> input)
{
	double sd = 0;
	unsigned tmpi = average(input);
	unsigned n = length(input);
	for (unsigned i=0;i<n;i++)
	{
		sd += (length(input[i])-tmpi)*(length(input[i])-tmpi);
	}
	sd /= n;
	sd = sqrt(sd);
	
	
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
int devide(T& input,T& less,T& great)
{
	int tmpi = average(input);
	for (unsigned i=0;i<length(input);i++)
	{
		if(length(input[i])<tmpi){appendValue(less,input[i]);}
		else{appendValue(great,input[i]);}
	}
	
	return 0;
}

int main(int argc, char *argv[])
{
	argv[1] = "/Informatik/Development/example.sam";
	seqan::String<seqan::String<char>> save;
	seqan::String<seqan::String<char>> less;
	seqan::String<seqan::String<char>> great;
	
	saminput(save,argv[1]);

	devide(save,less,great);

	std::cout << "less:" << std::endl;
	for (unsigned i=0;i<length(less);i++)
	{
		std::cout << less[i] << std::endl;
	}
	std::cout << "great:" << std::endl;
	for (unsigned i=0;i<length(great);i++)
	{
		std::cout << great[i] << std::endl;
	}

	
	return 0;
}