/*
Name: PEMer_Lite_execute
Author: Lars Zerbe <larszerbe@live.de>
Author: Stephan Peter <s.peter@fu-berlin.de>
License: GPL v3
Copyright: 2008-2013, FU Berlin
Status: under development
*/

#include "PEMer_Lite.h"

int main(int argc, char const ** argv)
{
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

	long double startTime = seqan::sysTime(), overallTime = startTime;

    // SAM-Datei einlesen
	saminput(save,header,options,stats);

    std::cout << "SAM input loading time: " << (seqan::sysTime() - startTime) << "s" << std::endl;

    startTime = seqan::sysTime();

    // Insertionen und Deletionen suchen
	devide(save,result,options,stats);

	std::cout << "indel searching time: " << (seqan::sysTime() - startTime) << "s" << std::endl;

    startTime = seqan::sysTime();

    // Clustering
    std::vector<std::vector<unsigned> > cluster;

	findCluster(result,cluster,candidate::insertion,stats);
	findCluster(result,cluster,candidate::deletion,stats);

    std::cout << "clustering time: " << (seqan::sysTime() - startTime) << "s" << std::endl;

    startTime = seqan::sysTime();

	std::string test = toCString(options.inputFile);
	test.erase(test.end()-4,test.end());

    // Ergebnisse in TSV-Datei schreiben
	output(cluster,header,result,toCString((seqan::CharString)(test+"_output.tsv")),stats);

    std::cout << "output writing time: " << (seqan::sysTime() - startTime) << "s" << std::endl;
    std::cout << "overall time: " << (seqan::sysTime() - overallTime) << "s" << std::endl;

	return 0;
}
