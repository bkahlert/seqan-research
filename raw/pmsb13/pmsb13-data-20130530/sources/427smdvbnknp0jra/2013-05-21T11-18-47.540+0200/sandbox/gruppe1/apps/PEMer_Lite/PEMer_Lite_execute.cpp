/*
Name: PEMer_Lite_execute
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
	seqan::String<seqan::CharString> header;
	saminput(save,header,options,stats);

	devide(save,result,options,stats);

    // Clustering
    std::vector<std::vector<unsigned> > cluster;

	findCluster(result,cluster,candidate::insertion,options,stats);
	findCluster(result,cluster,candidate::deletion,options,stats);

	std::string test = toCString(options.inputFile);
	test.erase(end(test)-4,end(test));

	output(cluster,header,result,toCString((seqan::CharString)(test+"_output.tsv")));

	return 0;
}
