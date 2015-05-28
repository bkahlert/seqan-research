#include <algorithm>
#include <iterator>
#include <iostream>
#include <seqan/bam_io.h>
#include <math.h>
#include <vector>
#include <fstream>


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


int main(int argc,char argv)
{	
	/*std::ofstream myfile;
	myfile.open ("/Informatik/Development/example.sam");
	myfile << "@HD	VN:1.3	SO:coordinate\n@SQ	SN:ref	LN:45\n@SQ	SN:ref2	LN:40\nr001	163	ref	7	30	8M4I4M1D3M	=	37	39	TTAGATAAAGAGGATACTG	*	XX:B:S,12561,2,20,112\nr002	0	ref	9	30	1S2I6M1P1I1P1I4M2I	*	0	0	AAAAGATAAGGGATAAA	*\nr003	0	ref	9	30	5H6M	*	0	0	AGCTAA	*\nr004	0	ref	16	30	6M14N1I5M	*	0	0	ATAGCTCTCAGC	*\nr003	16	ref	29	30	6H5M	*	0	0	TAGGC	*\nr001	83	ref	37	30	9M	=	7	-39	CAGCGCCAT	*";
	myfile.close();
	*/

	char *file = "/Informatik/Development/ecoli_indels_and_snps_reads.sam";

	candidate result,save;

	// BamStream open input stream
    seqan::BamStream SamInput(file);
    if (!isGood(SamInput))
    {
		std::cerr << "-ERROR- Could not open" << std::endl;
		return 1;
    }
	// record represent a record from a SAM-file
	seqan::BamAlignmentRecord record;

	std::cout << SamInput.bamIOContext._nameStore << std::endl;
	/*
	while (!atEnd(SamInput))
	{
		readRecord(record, SamInput);
		
		/*save.len.push_back(length(record.seq));
		save.ref.push_back(record.rID);
		save.pos.push_back(record.beginPos);
		save.type.push_back(candidate::u);
	}*/
	
	std::cout << "check" << std::endl;

    return 0;
}
