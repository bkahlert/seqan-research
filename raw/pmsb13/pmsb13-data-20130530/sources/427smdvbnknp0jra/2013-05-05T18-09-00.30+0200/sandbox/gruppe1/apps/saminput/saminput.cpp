/*
Name: saminput
Author: Lars Zerbe <larszerbe@live.de>
Maintainer: Lars Zerbe <larszerbe@live.de>
License: GPL v3
Copyright: 2008-2012, FU Berlin
Status: under development
*/

#include <iostream>
#include <seqan/bam_io.h>
#include <math.h>
#include <vector>

struct candidate
{
	enum form{u,d,i};	// form definiert den Type der Kandidaten ob u:undefiniert(0); d:deletion(1); i:insertion(2)

	std::vector<unsigned> pos;	// speichert die Positon
	std::vector<unsigned> ref;	// speichert die Referenz
	std::vector<unsigned> len;	// speichert die Laenge
	std::vector<form> type;		// speichert den Type
	unsigned av;				// speichert den Durchschnitt

	unsigned length()	// gibt die larnge aus *.length()
	{
		return len.size();
	}

	void clear()	// leert den Kandidaten
	{
		len.clear();
		pos.clear();
		ref.clear();
		type.clear();
	}
};

template <typename T>
T median(T& input)	// erstellt den median muss noch an struct candidate angepasst werden
{	
	std::nth_element( begin(input), begin(input) + length(input) / 2, end(input) ); 
	return *( begin(input) + length(input) / 2 ); 
}

unsigned average(candidate input)	// erzeugt den Durchschnitt und speichert ihm Objekt
{	
	unsigned tmpi=0;
	unsigned n = input.length();
	for (unsigned i=0;i<n;i++)
	{
		tmpi += input.len[i];
	}
	return tmpi/n;
}

unsigned analyse(candidate input)	//erzeugt die Standartabweichung
{

	double sd = 0;
	input.av=average(input);
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

int saminput(candidate &save, char *file)	// impotiert eine SAM-Datei und speichert sie im Format candidate
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
		save.type.push_back(candidate::u);
	}

    return 0;
}

int devide(candidate input,candidate &result)	// sucht nach insertions und deletions
{

	result.clear();

	unsigned sd = analyse(input);
	
	for (unsigned i=0;i<input.length();i++)
	{
		if(input.len[i]<input.av-sd)
		{
			
			result.len.push_back(input.len[i]);
			result.pos.push_back(input.pos[i]);
			result.ref.push_back(input.ref[i]);
			result.type.push_back(candidate::d);
		}
		if(input.len[i]>input.av+sd)
		{
			result.len.push_back(input.len[i]);
			result.pos.push_back(input.pos[i]);
			result.ref.push_back(input.ref[i]);
			result.type.push_back(candidate::i);
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
	// es werden alle typen an die Konsole ausgegen u=0, d=1, i=2
	if(unsigned n=result.length()!=0){
		std::cout << "candidates:" << std::endl;
		for (unsigned i=0;i<n;i++)
		{
			std::cout << result.type[i] << std::endl;
		}
	}
	
	return 0;
}