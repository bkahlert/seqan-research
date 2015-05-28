/*
Name: saminput
Author: Lars Zerbe <larszerbe@live.de>
Author: Stephan Peter <s.peter@fu-berlin.de>
Maintainer: Lars Zerbe <larszerbe@live.de>
License: GPL v3
Copyright: 2008-2012, FU Berlin
Status: under development
*/
#include <algorithm>
#include <iterator>
#include <iostream>
#include <seqan/bam_io.h>
#include <seqan/file.h>
#include <seqan/arg_parse.h>
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
	unsigned sd;				// speichert die Standardabweichung

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

struct modifier
{
    unsigned x,y;
    bool p;
    seqan::CharString inputFile;
    modifier() :
        x(1), y(1),p(false)
    {}
};

int median(candidate &input)	// erstellt den median muss noch an struct candidate angepasst werden
{
	std::nth_element( begin(input.len), begin(input.len) + input.length()/2,end(input.len));
	input.av=*(begin(input.len) + input.length()/2 );

	return 0;
}

unsigned average(candidate &input)	// erzeugt den Durchschnitt und speichert ihm Objekt
{
	unsigned tmpi=0;
	unsigned n = input.length();
	for (unsigned i=0;i<n;i++)
	{
		tmpi += input.len[i];
	}
	return tmpi/n;
}

seqan::ArgumentParser::ParseResult commands(modifier &options, int argc, char const ** argv)
{
    seqan::ArgumentParser parser("PEMer_Lite");
	// Optionen hinzufügen
    addOption(parser, seqan::ArgParseOption("i","input-file","Input SAM-file.",seqan::ArgParseArgument::INPUTFILE));
    setValidValues(parser, "input-file", "sam");
    setRequired(parser, "input-file");

    addOption(parser,seqan::ArgParseOption("x","SD-variable", "Variable to use for multiples of standard degression.",seqan::ArgParseArgument::INTEGER, "INT"));
	setMinValue(parser, "SD-variable", "1");
	setDefaultValue(parser,"SD-variable","1");

	addOption(parser,seqan::ArgParseOption("y","Del_SD-variable", "Variable to spezify multiples of standard degression for deletions.",seqan::ArgParseArgument::INTEGER, "INT"));
	setMinValue(parser, "Del_SD-variable", "1");
	setDefaultValue(parser,"Del_SD-variable","1");

	addOption(parser, seqan::ArgParseOption("p","Print", "Select to print temporary solutions."));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;
    // Extract option values.
	getOptionValue(options.x, parser, "SD-variable");
	getOptionValue(options.y, parser, "Del_SD-variable");
	options.p = isSet(parser, "Print");
	getOptionValue(options.inputFile, parser, "input-file");
    // If both to-uppercase and to-lowercase were selected then this is an error.
    
    return seqan::ArgumentParser::PARSE_OK;
}

int analyse(candidate &input)	//erzeugt die Standartabweichung
{

	median(input);
	input.sd=0;
	unsigned n = input.length();
	for (unsigned i=0;i<n;i++)
	{
		input.sd += (input.len[i]-input.av)*(input.len[i]-input.av);
	}
	input.sd /= n;
	input.sd = sqrt((double)input.sd);

	return 0;

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

int devide(candidate &input,candidate &result,modifier &options)	// sucht nach insertions und deletions
{
	result.clear();

	unsigned sd = analyse(input);

	for (unsigned i=0;i<input.length();i++)
	{
		if(input.len[i]<input.av-(options.x * sd))
		{

			result.len.push_back(input.len[i]);
			result.pos.push_back(input.pos[i]);
			result.ref.push_back(input.ref[i]);
			result.type.push_back(candidate::d);
		}
		if(input.len[i]>input.av+(options.y * sd))
		{
			result.len.push_back(input.len[i]);
			result.pos.push_back(input.pos[i]);
			result.ref.push_back(input.ref[i]);
			result.type.push_back(candidate::i);
		}
	}

	return 0;
}

int tsv(std::vector<std::vector<unsigned>> &input,candidate &ref, char *out)
 {	
	std::fstream outFile(out, std::ios::out);
	int k=0;

	if (!outFile) 
	{
		std::cerr << "-ERROR- Could not open output file" << std::endl;
		return 1;
	}
	else
	{	
		if(!input.empty())
		{
			for (unsigned i=0;i<input.size();i++)
			{	
				if(!input[i].empty())
				{
					outFile << "Cluster " << i-k+1 << ":\n";
					for(auto j=input[i].begin();j != input[i].end();j++)
					{
						outFile << "	" << ref.ref[*j] << "	" << ref.pos[*j] << "	" << ref.len[*j] << "	" << (candidate::form)ref.type[*j] << std::endl;
					}
				}
				else
				{
					k++;
				}
			}
		}
	}
	
	outFile.close();

	return 0;
 }

int findCluster(candidate &input, std::vector<std::vector<unsigned> > &dest, candidate::form indel)
{
    std::vector<unsigned> currentCluster;
    unsigned currentClusterEndPos;

    // Iterieren über alle Funde
    for(unsigned i = 0; i < input.length(); ++i)
    {
        // Insertionen bzw. Deletionen herausfiltern
        if(input.type[i] == indel)
        {
            if(currentCluster.empty())
            {
                // ersten Fund zu bisher leerem ersten Cluster hinzufügen
                currentCluster.push_back(i);
                currentClusterEndPos = input.pos[i] + input.len[i] - 1;
            }

            else
            {
                // Fall 1: Fund liegt innerhalb der Grenzen eines vorherigen Fundes
                if(input.pos[i] + input.len[i] - 1 <= currentClusterEndPos){
                    currentCluster.push_back(i);
				}

                else
                {
                    unsigned overlapLen = std::max((int)(currentClusterEndPos - input.pos[i] + 1), 0);

                    // Fall 2: Fund überlappt ausreichend mit dem Cluster
                    if(overlapLen / std::min(currentClusterEndPos - input.pos[currentCluster.back()] + 1, input.pos[i] + input.len[i] - 1) > 0.5){
                        currentCluster.push_back(i);
					}
                    // Fall 3: Fund überlappt nicht ausreichend mit dem Cluster
                    else
                    {
                        dest.push_back(currentCluster);
                        currentCluster.clear();
                        currentCluster.push_back(i);
                    }

                    currentClusterEndPos = input.pos[i] + input.len[i] - 1;
                }
            }
        }
    }

    dest.push_back(currentCluster);

	return 0;
}
 
int main(int argc, char const ** argv)
{	
	//char *file = "/Informatik/Development/test.sam";
	//argv[1] = "/Informatik/Development/example.sam";

	modifier options;
    seqan::ArgumentParser::ParseResult res = commands(options, argc, argv);
    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
	{
        return res == seqan::ArgumentParser::PARSE_ERROR;
	}

	candidate result,save;
	
	saminput(save,toCString(options.inputFile));
	
	devide(save,result,options);
	// es werden alle typen an die Konsole ausgegen u=0, d=1, i=2
	
		std::cout << "candidates:" << std::endl;
		for (unsigned i=0;i<result.length();i++)
		{
			std::cout << result.type[i] << "	" << result.pos[i] << "	" << result.ref[i] << std::endl;
		}
		
    // Clustering
    std::vector<std::vector<unsigned> > cluster;

    findCluster(result,cluster,candidate::i);
	findCluster(result,cluster,candidate::d);
	
	tsv(cluster,result,"/Informatik/Development/output.tsv");
	/*
		std::cout << "cluster2:" << std::endl;
		for (unsigned i=0;i<cluster.size();i++)
		{
		std::copy(cluster[i].begin(), cluster[i].end(), std::ostream_iterator<unsigned>(std::cout, " "));
		std::cout << std::endl;
		}
	*/

	return 0;
}
