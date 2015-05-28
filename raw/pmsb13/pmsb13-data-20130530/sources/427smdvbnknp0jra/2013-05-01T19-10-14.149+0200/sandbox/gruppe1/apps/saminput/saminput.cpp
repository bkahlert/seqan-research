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
		}
		return av/len.size();
	}
};

template <typename T>
T median(T& input)
{	
	std::nth_element( begin(input), begin(input) + length(input) / 2, end(input) ); 
	return *( begin(input) + length(input) / 2 ); 
}

template <typename T>
int average(candidate input)
{	
	int tmpi=0;
	for (unsigned i=0;i<input.numOfCan();i++)
	{
		tmpi += length(input[i]);
	}
	return tmpi/length(input);
}

unsigned analyse(candidate input)
{

	double sd = 0;
	input.average();
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
	unsigned sd = analyse(input);
	
	for (unsigned i=0;i<input.length();i++)
	{
		if((input.len[i]<sd-input.av) || (input.len[i]>sd+input.av))
		{
			/*result.len[i]=input.len[i];
			result.pos[i]=input.pos[i];
			result.ref[i]=input.ref[i];
		*/}
	}
	
	return 0;
}

int main(int argc, char *argv[])
{
	argv[1] = "/Informatik/Development/example.sam";
	
	candidate result,save;
		
	saminput(save,argv[1]);
	
	devide(save,result);
	/*
	if(unsigned n=result.length()!=0){
		std::cout << "candidates:" << std::endl;
		for (unsigned i=0;i<n;i++)
		{
			std::cout << result.ref[i] << std::endl;
		}
	}
	*/
	return 0;
}