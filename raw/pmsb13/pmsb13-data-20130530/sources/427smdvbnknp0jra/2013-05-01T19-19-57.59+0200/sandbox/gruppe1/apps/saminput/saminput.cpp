#include <iostream>
#include <seqan/bam_io.h>
#include <math.h>
#include <vector>

struct candidate
{
	enum form{d,i};

	std::vector<unsigned> pos;
	std::vector<unsigned> ref;
	std::vector<unsigned> len;
	unsigned av;
	form type;

	unsigned length()
	{
		return len.size();
	}

	unsigned average()
	{
		av = 0;
		for (unsigned i=0;i<len.size();i++)
		{
			av += len[i];
			std::cout << av << std::endl;
		}
		return (av/len.size());
	}
	void clear()
	{
		len.clear();
		pos.clear();
		ref.clear();
	}
};

template <typename T>
T median(T& input)
{	
	std::nth_element( begin(input), begin(input) + length(input) / 2, end(input) ); 
	return *( begin(input) + length(input) / 2 ); 
}


unsigned average(candidate input)
{	
	unsigned tmpi=0;
	unsigned n = input.length();
	for (unsigned i=0;i<n;i++)
	{
		tmpi += input.len[i];
	}
	return tmpi/n;
}

unsigned analyse(candidate input)
{

	double sd = 0;
	input.av=average(input);
	std::cout << input.av << std::endl;
	unsigned n = input.length();
	for (unsigned i=0;i<n;i++)
	{
		sd += (input.len[i]-input.av)*(input.len[i]-input.av);
	}
	sd /= n;
	sd = sqrt(sd);

	(unsigned)sd;
	
	return (unsigned)sd;
	
}

int saminput(candidate &save, char *file)
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
		save.len.push_back(length(record.seq));
		save.ref.push_back(record.rID);
		save.pos.push_back(record.beginPos);
	}

    return 0;
}


int devide(candidate input,candidate result)
{

	result.clear();

	unsigned sd = analyse(input);
	std::cout << sd << std::endl;
	for (unsigned i=0;i<input.length();i++)
	{
		if((input.len[i]<sd-input.av) || (input.len[i]>sd+input.av))
		{
			result.len.push_back(input.len[i]);
			result.pos.push_back(input.pos[i]);
			result.ref.push_back(input.ref[i]);
		}
	}
	
	return 0;
}

int main(int argc, char *argv[])
{
	argv[1] = "/Informatik/Development/example.sam";
	
	candidate result,save;
		
	saminput(save,argv[1]);
	
	devide(save,result);
	
	if(unsigned n=result.length()!=0){
		std::cout << "candidates:" << std::endl;
		for (unsigned i=0;i<n;i++)
		{
			std::cout << result.ref[i] << std::endl;
		}
	}
	
	return 0;
}