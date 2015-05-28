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

	addOption(parser,seqan::ArgParseOption("s","standard_degression", "Variable to use for standard degression for length of fragments.",seqan::ArgParseArgument::INTEGER, "INT"));
	setMinValue(parser, "standard_degression", "1");
	
	addOption(parser,seqan::ArgParseOption("e","expected_value", "Variable to use for the expected value for length of fragments.",seqan::ArgParseArgument::INTEGER, "INT"));
	setMinValue(parser, "expected_value", "1");
		
	addOption(parser, seqan::ArgParseOption("p","Print", "Select to print temporary solutions."));

    // verarbeitet die Kommandozeile
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Überprüfung der richtigen Ausführung von commands()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // extrahiert Werte aus der parser.
	getOptionValue(options.x, parser, "SD-variable");
	getOptionValue(options.y, parser, "Del_SD-variable");
	getOptionValue(options.e, parser, "expected_value");
	getOptionValue(options.sd, parser, "standard_degression");
	options.p = isSet(parser, "Print");
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

int saminput(candidates &save, modifier &options, statisticals &stats)	// impotiert eine SAM-Datei und speichert sie im Format candidates
{	
	candidate tmp;
	medianTree median;
	
    // BamStream öffnet Eingangsdatenstrom
    seqan::BamStream SamInput(toCString(options.inputFile));
    if (!isGood(SamInput))
    {
		std::cerr << "-ERROR- Could not open" << std::endl;
		return 1;
    }
	/*
	for (unsigned i = 0; i < length(SamInput.header.sequenceInfos); ++i)
	{
		header << SamInput.header.sequenceInfos[i].i1 << "\t" << SamInput.header.sequenceInfos[i].i2 << "\n";
	}*/
	// record enthält die inhalte aus einer SAM-Datei
	seqan::BamAlignmentRecord record;

	while (!atEnd(SamInput))
	{
		readRecord(record, SamInput);
		if(record.tLen>0)
		{
			median.add(record.tLen);
			tmp.addvalues(record.rID,record.beginPos,record.tLen,candidate::undefine);
			save.push_back(tmp);
		}
	}
	if(options.e==NULL)
	{
		stats.e_value=median.getMedian();
	}
	else
	{
		stats.e_value=options.e;
	}
	if(options.sd==NULL)
	{
		stats.e_value=generateStandardDegression(save,stats.e_value);
	}
	else
	{
		stats.standardDegression=options.sd;
	}
	
    return 0;
}

int devide(candidates &input, candidates &result, modifier &options, statisticals &stats)	// sucht nach insertions und deletions
{
	for (unsigned i=0;i<input.size();i++)
	{
		if(input[i].len<stats.e_value-(options.x * stats.standardDegression))
		{

			result.push_back(input[i]);
			result[i].type=candidate::deletion;
		}
		if(input[i].len>stats.e_value+(options.y * stats.standardDegression))
		{
			result.push_back(input[i]);
			result[i].type=candidate::insertion;
		}
	}

	return 0;
}

int output(std::vector<std::vector<unsigned>> &input,candidates &ref, char *out)
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
			for(unsigned j=0;j<input.size();j++)
			{
				if(!input[j].empty())
				{
					unsigned targetSize=0;
					unsigned maxBegin=ref[0].pos;
					unsigned minEnd=ref[0].pos+ref[0].len;
					unsigned n=input[j].size();
					outFile << "Cluster: " << std::endl;

					for(auto k=input[j].begin();k != input[j].end();k++)
					{
						targetSize += ref[*k].len;
						
						if(ref[*k].pos > maxBegin){maxBegin=ref[*k].pos;}
						if((ref[*k].pos+ref[*k].len) < minEnd){minEnd=(ref[*k].pos+ref[*k].len);}
					}
					targetSize /= n;
					outFile << ref[0].ref << "\t" << (maxBegin+minEnd)/2<< "\t" << targetSize << "\t" << maxBegin << "\t" << minEnd << "\t" << (candidate::form)ref[0].type << std::endl;
				}
			}
		}
	}

	outFile.close();

	return 0;
}
 
int findCluster(candidates &input, std::vector<std::vector<unsigned> > &dest, candidate::form indel)
{
    std::vector<unsigned> currentCluster;
    unsigned currentClusterEndPos;

    // Iterieren über alle Funde
	for(unsigned i = 0; i < input.size(); ++i)
    {
        // Insertionen bzw. Deletionen herausfiltern
        if(input[i].type == indel)
        {
            if(currentCluster.empty())
            {
                // ersten Fund zu bisher leerem ersten Cluster hinzufügen
                currentCluster.push_back(i);
                currentClusterEndPos = input[i].pos + input[i].len - 1;
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

                    // Fall 2: Fund überlappt ausreichend mit dem Cluster
                    if(overlapLen / std::min(currentClusterEndPos - input[currentCluster.back()].pos + 1, input[i].pos + input[i].len - 1) > 0.5){
                        currentCluster.push_back(i);
					}
                    // Fall 3: Fund überlappt nicht ausreichend mit dem Cluster
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

    dest.push_back(currentCluster);

	return 0;
}
 
