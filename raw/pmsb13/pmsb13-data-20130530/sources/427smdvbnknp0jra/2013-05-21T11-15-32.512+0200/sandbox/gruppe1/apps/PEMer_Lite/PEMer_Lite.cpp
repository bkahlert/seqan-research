/*
Name: PEMer_Lite
Author: Lars Zerbe <larszerbe@live.de>
Author: Stephan Peter <s.peter@fu-berlin.de>
Maintainer: Lars Zerbe <larszerbe@live.de>
License: GPL v3
Copyright: 2008-2012, FU Berlin
Status: under development
*/

#include "PEMer_Lite.h"


bool candidateSortFunction(candidate x, candidate y)
{
	try{return (x.pos < y.pos);}
	catch(std::exception allgemeineException){return false;}
}

int samSort(candidates &save)
{
	std::sort(save.begin(),save.end(),candidateSortFunction);

	return 0;
}

double coverage(unsigned l,unsigned n, unsigned g)
{
	return (1- std::exp((double)(l*n/g)));
	/*
	G is the  genome length
	L is the read length
	N is the number of reads
	*/
}

seqan::ArgumentParser::ParseResult commands(modifier &options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("PEMer_Lite");

	// Optionen hinzuf�gen
    addOption(parser, seqan::ArgParseOption("f","input-file","Input SAM-file.",seqan::ArgParseArgument::INPUTFILE));
    setValidValues(parser, "input-file", "sam");
    setRequired(parser, "input-file");

    addOption(parser,seqan::ArgParseOption("x","SD-variable", "Variable to use for multiples of standard degression.",seqan::ArgParseArgument::INTEGER, "INT"));
	setMinValue(parser, "SD-variable", "1");
	setDefaultValue(parser,"SD-variable","1");

	addOption(parser,seqan::ArgParseOption("y","Del_SD-variable", "Variable to spezify multiples of standard degression for deletions.",seqan::ArgParseArgument::INTEGER, "INT"));
	setMinValue(parser, "Del_SD-variable", "1");
	setDefaultValue(parser,"Del_SD-variable","1");

	addOption(parser,seqan::ArgParseOption("i","max_insertion_size", "This value spezify the maximum insertion size.",seqan::ArgParseArgument::INTEGER, "INT"));
	setMinValue(parser, "max_insertion_size", "1");
	setDefaultValue(parser,"max_insertion_size","100000");

	addOption(parser,seqan::ArgParseOption("d","max_deletion_size", "This value spezify the maximum deletion size.",seqan::ArgParseArgument::INTEGER, "INT"));
	setMinValue(parser, "max_deletion_size", "1");
	setDefaultValue(parser,"max_deletion_size","1000000");

	addOption(parser,seqan::ArgParseOption("c","coverage", "Variable to use for p value coverage.",seqan::ArgParseArgument::INTEGER, "INT"));
	setMinValue(parser, "coverage", "0");


	addOption(parser,seqan::ArgParseOption("sd","standard_degression", "Variable to use for standard degression for length of fragments.",seqan::ArgParseArgument::INTEGER, "INT"));
	setMinValue(parser, "standard_degression", "1");

	addOption(parser,seqan::ArgParseOption("e","expected_value", "Variable to use for the expected value for length of fragments.",seqan::ArgParseArgument::INTEGER, "INT"));
	setMinValue(parser, "expected_value", "1");

	addOption(parser, seqan::ArgParseOption("p","Print", "Select to print temporary solutions."));

	addOption(parser, seqan::ArgParseOption("s","Sort_Position", "Select to sort for Position."));

    // verarbeitet die Kommandozeile
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // �berpr�fung der richtigen Ausf�hrung von commands()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // extrahiert Werte aus der parser.
	getOptionValue(options.x, parser, "SD-variable");
	getOptionValue(options.y, parser, "Del_SD-variable");
	getOptionValue(options.e, parser, "expected_value");
	getOptionValue(options.sd, parser, "standard_degression");
	getOptionValue(options.i, parser, "max_insertion_size");
	getOptionValue(options.d, parser, "max_deletion_size");
	getOptionValue(options.c, parser, "coverage");
	options.isC = isSet(parser, "coverage");
	options.p = isSet(parser, "Print");
	options.s = isSet(parser, "Sort_Position");
	getOptionValue(options.inputFile, parser, "input-file");

    return seqan::ArgumentParser::PARSE_OK;
}

unsigned generateStandardDegression(candidates &input,int e_value)	//erzeugt die Standartabweichung
{
	unsigned standardDegression=0;

	unsigned n = input.size();
	for (unsigned i=0;i<n;i++)
	{
		standardDegression += (input[i].len-e_value)*(input[i].len-e_value);
	}
	standardDegression /= n;
	standardDegression = (unsigned)sqrt((double)standardDegression);

	return standardDegression;
}

int saminput(candidates &save, seqan::String<seqan::CharString> &header, modifier &options, statisticals &stats)	// impotiert eine SAM-Datei und speichert sie im Format candidates
{
	candidate tmp;
	medianTree median;
	unsigned n=0;
	
    // BamStream �ffnet Eingangsdatenstrom
    seqan::BamStream SamInput(toCString(options.inputFile));
    if (!isGood(SamInput))
    {
		std::cerr << "-ERROR- Could not open" << std::endl;
		return 1;
    }

	for(unsigned i=0; i<length(SamInput.header.sequenceInfos);i++)
	{
		seqan::appendValue(header,SamInput.header.sequenceInfos[i].i1)
	}
	// record enth�lt die inhalte aus einer SAM-Datei
	seqan::BamAlignmentRecord record;

	while (!atEnd(SamInput))
	{
		readRecord(record, SamInput);
		n++;
		if(record.tLen>0)
		{
			median.add(record.tLen);
			tmp.addvalues(record.rID,record.beginPos,record.tLen,candidate::undefine);
			save.push_back(tmp);
		}
	}
	
	if(options.e==0)
	{
		stats.e_value=median.getMedian();
	}
	else
	{
		stats.e_value=options.e;
	}
	if(options.sd==0)
	{
		stats.standardDegression=generateStandardDegression(save,stats.e_value);
	}
	else
	{
		stats.standardDegression=options.sd;
	}
	stats.coverage = coverage(length(record.seq),n,SamInput.header.sequenceInfos[0].i2);
	
    return 0;
}

int devide(candidates &input, candidates &result, seqan::String<seqan::CharString> &header, modifier &options, statisticals &stats)	// sucht nach insertions und deletions
{
	for (unsigned i=0;i<input.size();i++)
	{
		if(input[i].len<(stats.e_value-(options.x * stats.standardDegression)))
		{
			input[i].type=candidate::deletion;
			result.push_back(input[i]);
		}else
			if(input[i].len>(stats.e_value+(options.y * stats.standardDegression)))
			{
				input[i].type=candidate::insertion;
				result.push_back(input[i]);
			}
	}
	if(options.s){samSort(result);}

	return 0;
}

int output(std::vector<std::vector<unsigned> > &input,candidates &ref, char *out)
{
	std::fstream outFile(out, std::ios::out);

	if (!outFile)
	{
		std::cerr << "-ERROR- Could not open output file" << std::endl;
		return 1;
	}
	else
	{
		outFile << "ID\tPOS\tLEN\tBegin\tEND\tTYPE " << std::endl;
		if(!input.empty())
		{
			for (unsigned i=0;i<input.size();i++)
			{
				unsigned targetSize=0;
				unsigned maxBegin=0;
				unsigned minEnd=0;
				unsigned n=input[i].size();

				if(!input[i].empty())
				{
					for(std::vector<unsigned>::iterator j=input[i].begin();j != input[i].end();j++)
					{
						targetSize += ref[*j].len;

						if(ref[*j].pos > maxBegin){maxBegin=ref[*j].pos;}
						if((ref[*j].pos+ref[*j].len) < minEnd||minEnd==0){minEnd=(ref[*j].pos+ref[*j].len);}
					}
				}
				if(n!=0)
				{
					targetSize /= n;

					outFile << header[ref[*input[i].begin()].ref] << "\t" << (maxBegin+minEnd)/2<< "\t" << targetSize << "\t" << maxBegin << "\t" << minEnd << "\t" << (candidate::form)ref[*input[i].begin()].type << std::endl;
				}
			}
		}
	}

	outFile.close();

	return 0;
}

void findCluster(candidates &input, std::vector<std::vector<unsigned> > &dest, candidate::form indel, modifier &options, statisticals &stats)
{
    std::vector<unsigned> currentCluster;
    unsigned currentClusterEndPos;
    bool emptyFlag = true;
	double coverage;
	if(options.isC){coverage=(double)options.c;}else{coverage=stats.coverage;}
		
	// Iterieren �ber alle Funde
	for(unsigned i = 0; i < input.size(); ++i)
    {
        // Insertionen bzw. Deletionen herausfiltern
        if(input[i].type == indel)
        {
            if(currentCluster.empty())
            {
                // ersten Fund zu bisher leerem ersten Cluster hinzuf�gen
                currentCluster.push_back(i);
                currentClusterEndPos = input[i].pos + input[i].len - 1;
                emptyFlag = false;
            }

            else
            {
                // Fall 1: Fund liegt innerhalb der Grenzen eines vorherigen Fundes
                if(input[i].pos + input[i].len - 1 <= currentClusterEndPos){
                    currentCluster.push_back(i);
				}

                else
                {
                    unsigned overlapLen = std::max((int)(currentClusterEndPos - input[i].pos + 1), 0);

                    // Fall 2: Fund �berlappt ausreichend mit dem Cluster
                    if((double)overlapLen / std::min(currentClusterEndPos - input[currentCluster.back()].pos + 1, input[i].pos + input[i].len - 1) > coverage){
                        currentCluster.push_back(i);
					}

                    // Fall 3: Fund �berlappt nicht ausreichend mit dem Cluster
                    else
                    {
                        dest.push_back(currentCluster);
                        currentCluster.clear();
                        currentCluster.push_back(i);
                    }

                    currentClusterEndPos = input[i].pos + input[i].len - 1;
                }
            }
        }
    }

    if(!emptyFlag)
        dest.push_back(currentCluster);
}

/*
void findCluster(candidates &input, std::vector<std::vector<unsigned> > &dest, candidate::form indel, double coverage)
{
    std::vector<unsigned> currentCluster;
    unsigned currentClusterEndPos;
    bool emptyFlag = true;

    // Iterieren �ber alle Funde
	for(unsigned i = 0; i < input.size(); ++i)
    {
        // Insertionen bzw. Deletionen herausfiltern
        if(input[i].type == indel)
        {
            if(currentCluster.empty())
            {
                // ersten Fund zu bisher leerem ersten Cluster hinzuf�gen
                currentCluster.push_back(i);
                currentClusterEndPos = input[i].pos + input[i].len - 1;
                emptyFlag = false;
            }

            else
            {
                // Fall 1: Fund liegt innerhalb der Grenzen eines vorherigen Fundes
                if(input[i].pos + input[i].len - 1 <= currentClusterEndPos){
                    currentCluster.push_back(i);
				}

                else
                {
                    unsigned overlapLen = std::max((int)(currentClusterEndPos - input[i].pos + 1), 0);

                    // Fall 2: Fund �berlappt ausreichend mit dem Cluster
                    if((double)overlapLen / std::min(currentClusterEndPos - input[currentCluster.back()].pos + 1, input[i].pos + input[i].len - 1) > coverage){
                        currentCluster.push_back(i);
					}

                    // Fall 3: Fund �berlappt nicht ausreichend mit dem Cluster
                    else
                    {
                        dest.push_back(currentCluster);
                        currentCluster.clear();
                        currentCluster.push_back(i);
                    }

                    currentClusterEndPos = input[i].pos + input[i].len - 1;
                }
            }
        }
    }

    if(!emptyFlag)
        dest.push_back(currentCluster);
}

*/