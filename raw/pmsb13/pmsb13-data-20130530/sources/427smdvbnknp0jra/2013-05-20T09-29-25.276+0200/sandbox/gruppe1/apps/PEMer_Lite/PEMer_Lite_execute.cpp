/*
Name: saminput
Author: Lars Zerbe <larszerbe@live.de>
Author: Stephan Peter <s.peter@fu-berlin.de>
Maintainer: Lars Zerbe <larszerbe@live.de>
License: GPL v3
Copyright: 2008-2012, FU Berlin
Status: under development
*/

#include "PEMer_Lite.h"

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

	candidates result,save;
	statisticals stats;
		
	saminput(save,options,stats);
	
	devide(save,result,options,stats);

	samSort(save);
	for(unsigned i=0;i<save.size();i++)
	{
		std::cout << save[i].len<<"\t" << save[i].pos << "\t"<<save[i].type;
	}
    // Clustering
    std::vector<std::vector<unsigned> > cluster;

	findCluster(result,cluster,candidate::insertion);
	findCluster(result,cluster,candidate::deletion);
	
	std::string test = toCString(options.inputFile);
	test.erase(end(test)-4,end(test));
	
	output(cluster,result,toCString((seqan::CharString)(test+"_output.tsv")));
	
	return 0;
}
