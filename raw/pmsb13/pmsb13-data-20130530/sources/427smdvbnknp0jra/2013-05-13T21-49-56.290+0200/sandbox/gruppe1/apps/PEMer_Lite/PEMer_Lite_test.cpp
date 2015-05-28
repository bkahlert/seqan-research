#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include "PEMer_Lite.h"



int median(candidate &input)	// erstellt den median muss noch an struct candidate angepasst werden
{
	std::nth_element( begin(input.len), begin(input.len) + input.length()/2,end(input.len));
	input.e=*(begin(input.len) + input.length()/2 );

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

int analyse(candidate &input)	//erzeugt die Standartabweichung
{

	median(input);
	input.sd=0;
	unsigned n = input.length();
	for (unsigned i=0;i<n;i++)
	{
		input.sd += (input.len[i]-input.e)*(input.len[i]-input.e);
	}
	input.sd /= n;
	input.sd = sqrt((double)input.sd);

	return 0;

}

int saminput(candidate &save, char *file)	// impotiert eine SAM-Datei und speichert sie im Format candidate
{
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
	if(options.sd==0)
	{
		analyse(input);
	}
	else
	{
		input.e=options.e;
		input.sd=options.sd;
	}
	for (unsigned i=0;i<input.length();i++)
	{
		if(input.len[i]<input.e-(options.x * input.sd))
		{

			result.len.push_back(input.len[i]);
			result.pos.push_back(input.pos[i]);
			result.ref.push_back(input.ref[i]);
			result.type.push_back(candidate::d);
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



SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_input)
{
	candidate test;
	char *file ="/Informatik/Development/test.sam";
	SEQAN_ASSERT(!saminput(test,file));
}

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_input2)
{
	candidate test;
	char *file =" ";
	SEQAN_ASSERT(saminput(test,file));
}

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_input3)
{
	candidate test;
	char *file ="2579816ß15362466183b513ü6ß148z5!§%%!&°%(!%§(%§c465cn4ü457346511425129562e85z148z5!§%%!&°%(!%§(%§&&&&&(";
	SEQAN_ASSERT(saminput(test,file));
}

SEQAN_DEFINE_TEST(test_my_app_PEMer_Lite_devide)
{
	candidate test,result;
	modifier opt;

	for(unsigned i=1;i<6;i++)
	{
		test.len.push_back(10);
		test.pos.push_back(i);
		test.ref.push_back(i);
	}
	test.len.push_back(3);
		test.pos.push_back(6);
		test.ref.push_back(6);
	test.len.push_back(9);
		test.pos.push_back(7);
		test.ref.push_back(7);
	test.len.push_back(11);
		test.pos.push_back(8);
		test.ref.push_back(8);
	test.len.push_back(1000);
		test.pos.push_back(9);
		test.ref.push_back(9);
	test.len.push_back(200);
		test.pos.push_back(10);
		test.ref.push_back(10);

	test.e=10;
	test.sd=1;
	opt.x=1;
	opt.y=1;

	devide(test,result,opt);

	SEQAN_ASSERT_EQ(result.type[0],1);
	SEQAN_ASSERT_EQ(result.type[1],1);
	SEQAN_ASSERT_EQ(result.type[2],1);
	
	
	SEQAN_ASSERT_EQ(result.ref[0],6);
	SEQAN_ASSERT_EQ(result.ref[1],9);
	SEQAN_ASSERT_EQ(result.ref[2],10);
}


SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
   SEQAN_CALL_TEST(test_my_app_PEMer_Lite_input);
   SEQAN_CALL_TEST(test_my_app_PEMer_Lite_input2);
   SEQAN_CALL_TEST(test_my_app_PEMer_Lite_input3);
   SEQAN_CALL_TEST(test_my_app_PEMer_Lite_devide);
}
SEQAN_END_TESTSUITE