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

/*
bool compareCandidateLength(candidate x,candidate y)
{
	return (x.len < y.len);
}

unsigned median(candidates &input)	// erstellt den median muss noch an struct candidate angepasst werden
{
	std::nth_element(input.begin(),input.begin() + input.size()/2,input.end(),compareCandidateLength);
	return input[input.size()/2].len;
	
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

statisticals generateStandardDegression(candidates &input)	//erzeugt die Standartabweichung
{
	statisticals output;

	output.e_value=median(input);

	unsigned n = input.size();
	for (unsigned i=0;i<n;i++)
	{
		output.standardDegression += (input[i].len-output.e_value)*(input[i].len-output.e_value);
	} 
	output.standardDegression /= n;
	output.standardDegression = sqrt((double)output.standardDegression);

	return output;
}

int saminput(candidates &save, char *file)	// impotiert eine SAM-Datei und speichert sie im Format candidate
{
	candidate tmp;

    // BamStream öffnet Eingangsdatenstrom
    seqan::BamStream SamInput(file);
    if (!isGood(SamInput))
    {
		std::cerr << "-ERROR- Could not open" << std::endl;
		return 1;
    }
	// record enthält die inhalte aus einer SAM-Datei
	seqan::BamAlignmentRecord record;

	while (!atEnd(SamInput))
	{
		readRecord(record, SamInput);
		tmp.addvalues(record.rID,record.beginPos,record.tLen,candidate::undefine);
		if(tmp.len>0){save.push_back(tmp);}
	}

    return 0;
}

int devide(candidates &input,candidates &result,modifier &options)	// sucht nach insertions und deletions
{
	statisticals stats;
	if(options.sd==0)
	{
		stats=generateStandardDegression(input);
	}
	else
	{
		stats.e_value=options.e;
		stats.standardDegression=options.sd;
	}
	for (unsigned i=0;i<input.size();i++)
	{
		if(input[i].len<stats.e_value-(options.x * stats.standardDegression))
		{

			result.push_back(input[i])
		}
		if(input.len[i]>input.e+(options.y * input.sd))
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
 */
int main(int argc, char const ** argv)
{	/*
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
	
	std::string test = toCString(options.inputFile);
	test.erase(end(test)-4,end(test));
	
	tsv(cluster,result,toCString((seqan::CharString)(test+"_output.tsv")));
	*/

	medianTree median;
	median.add(3);
	median.add(1);
	median.add(2);
	median.add(3);
	median.add(3);
	std::cout << median.getMedian();
	return 0;
}
