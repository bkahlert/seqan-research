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

seqan::String<seqan::String<unsigned>> analyse(seqan::String<seqan::String<char>> input)
{

	double sd = 0;
	unsigned tmpi = average(input);
	unsigned n = length(input);
	seqan::String<seqan::String<unsigned>> result;
	for (unsigned i=0;i<n;i++)
	{
		sd += (length(input[i])-tmpi)*(length(input[i])-tmpi);
	}
	sd /= n;
	sd = sqrt(sd);

	appendValue(result,(unsigned)sd);
	appendValue(result,tmpi);

	return result;
	
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
int devide(T& input,T& result,unsigned av)
{
	seqan::String<seqan::String<unsigned>> tmp = analyse(input);

	for (unsigned i=0;i<length(input);i++)
	{
		if(length(input[i])<tmp[2]-tmp[1] || length(input[i])<tmp[2]+tmp[1]){appendValue(result,input[i]);}
	}
	
	return 0;
}

int main(int argc, char *argv[])
{
	argv[1] = "/Informatik/Development/example.sam";
	seqan::String<seqan::String<char>> save;
	seqan::String<seqan::String<char>> result;
		
	saminput(save,argv[1]);

	
	if(length(result)!=0){
		std::cout << "candidates:" << std::endl;
		for (unsigned i=0;i<length(result);i++)
		{
			std::cout << result[i] << std::endl;
		}
		
	std::cout << analyse(save) << std::endl;
	return 0;
}