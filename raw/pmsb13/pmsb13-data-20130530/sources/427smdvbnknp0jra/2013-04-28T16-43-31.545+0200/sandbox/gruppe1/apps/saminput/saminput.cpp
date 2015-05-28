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

int analyse(seqan::String<seqan::String<char>> input,seqan::String<seqan::String<unsigned>> &result)
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

	appendValue(result,(unsigned)sd);
	appendValue(result,tmpi);

	return 0;
	
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
		appendValue(save,record.seq);		
	}

    return 0;
}

template <typename T>
int devide(T& input,T& result)
{
	seqan::String<seqan::String<unsigned>> tmp;
	//analyse(input,tmp);
	appndValue(tmp,5);
	appndValue(tmp,11);
	for (unsigned i=0;i<length(input);i++)
	{
		if(length(input[i])<tmp[2]-tmp[1] || length(input[i])>tmp[2]+tmp[1]){appendValue(result,input[i]);}
	}
	
	return 0;
}

int main(int argc, char *argv[])
{
	argv[1] = "/Informatik/Development/example.sam";
	seqan::String<seqan::String<char>> save;
	seqan::String<seqan::String<char>> result;
		
	saminput(save,argv[1]);
	
	//devide(save,result);
	for(int k=0;k<length(result);k++)
	{
		std::cout << result[k] << std::endl;
	}
	if(unsigned n=length(result)!=0){
		std::cout << "candidates:" << std::endl;
		for (unsigned i=0;i<n;i++)
		{
			std::cout << result[i] << std::endl;
		}
	}
	
	return 0;
}