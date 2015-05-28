/*
Name: PEMer_Lite
Author: Lars Zerbe <larszerbe@live.de>
Author: Stephan Peter <s.peter@fu-berlin.de>
License: GPL v3
Copyright: 2008-2013, FU Berlin
Status: under development
*/

#include "PEMer_Lite.h"

bool candidateSortFunction(candidate x, candidate y)
{	// Sortierverfahren den struct candidate.
	try{return (x.pos < y.pos);}
	catch(std::exception allgemeineException){return false;}
}

int samSort(candidates &save)
{	// Sortier aufruf f�r candidates, welche nach ihrer Position im Genome sortiert werden.
	std::sort(save.begin(),save.end(),candidateSortFunction);

	return 0;
}

int coverage(modifier &options, statisticals &stats, unsigned l,unsigned n, unsigned g)
{
	if(options.isC)
	{
		stats.fold_sequence_coverage=1/std::exp((double)options.c);
		stats.coverage=(double)options.c;
	}else
	{
		stats.fold_sequence_coverage=1/std::exp((double)(l*n/g));
		stats.coverage=(double)(l*n/g);
	}
	
	return 0;
	/*
	G L�nge des Genomes
	L L�nge der Reads
	N Anzahl der Reads
	Daraus wird die coverage und die "fold sequence coverage" berechnet, welche die Abdeckung darstellt in Werten zwischen 0 und 1.
	*/
}

seqan::ArgumentParser::ParseResult commands(modifier &options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("PEMer_Lite");

	// Es werden Optionen f�r die Kommandozeile hinzugef�gt.
    addOption(parser, seqan::ArgParseOption("f","input-file","Input SAM-file.",seqan::ArgParseArgument::INPUTFILE));
    setValidValues(parser, "input-file", "sam");
    setRequired(parser, "input-file");

    addOption(parser,seqan::ArgParseOption("x","SD-variable", "Variable to use for multiples of standard degression.",seqan::ArgParseArgument::INTEGER, "INT"));
	setMinValue(parser, "SD-variable", "1");
	setDefaultValue(parser,"SD-variable","2");

	addOption(parser,seqan::ArgParseOption("y","Del_SD-variable", "Variable to spezify multiples of standard degression for deletions.",seqan::ArgParseArgument::INTEGER, "INT"));
	setMinValue(parser, "Del_SD-variable", "1");
	setDefaultValue(parser,"Del_SD-variable","2");

	addOption(parser,seqan::ArgParseOption("c","coverage", "Variable to use for p value coverage.",seqan::ArgParseArgument::INTEGER, "INT"));
	setMinValue(parser, "coverage", "0");

	addOption(parser,seqan::ArgParseOption("sd","standard_degression", "Variable to use for standard degression for length of fragments.",seqan::ArgParseArgument::INTEGER, "INT"));
	setMinValue(parser, "standard_degression", "1");

	addOption(parser,seqan::ArgParseOption("e","expected_value", "Variable to use for the expected value for length of fragments.",seqan::ArgParseArgument::INTEGER, "INT"));
	setMinValue(parser, "expected_value", "1");

	addOption(parser, seqan::ArgParseOption("p","Print", "Select to print temporary solutions."));

	addOption(parser, seqan::ArgParseOption("s","Sort_Position", "Select to sort for Position."));

    // Die Kommandozeile wird verarbeitet.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // �berpr�fung der richtigen Ausf�hrung von commands()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extrahiert Werte aus der Kommandozeile.
	getOptionValue(options.x, parser, "SD-variable");
	getOptionValue(options.y, parser, "Del_SD-variable");
	getOptionValue(options.e, parser, "expected_value");
	getOptionValue(options.sd, parser, "standard_degression");
	getOptionValue(options.c, parser, "coverage");
	options.isC = isSet(parser, "coverage");
	options.p = isSet(parser, "Print");
	options.s = isSet(parser, "Sort_Position");
	getOptionValue(options.inputFile, parser, "input-file");

    return seqan::ArgumentParser::PARSE_OK;
}

unsigned generateStandardDegression(candidates &input,int e_value)	// Berechnet die Standartabweichung
{
	unsigned standardDegression=0;	// Die Standardabweichung wird mit 0 initialisiert.

	unsigned n = input.size();	// Es wird �ber alle Elemente der Eingabe iteriert.
	for (unsigned i=0;i<n;i++)
	{
		standardDegression += (input[i].len-e_value)*(input[i].len-e_value);	// Die Standardabweichung wird um (x-E(X))^2 erh�ht.
	}
	standardDegression /= n;	// Die Standardabweichung wird durch die Anzahl der Elemnte der Eingabe geteilt.
	standardDegression = (unsigned)sqrt((double)standardDegression);	// Es wird die Wurzel der Standardabweichung gebildet.

	return standardDegression;
}

int saminput(candidates &save, seqan::String<seqan::CharString> &header, modifier &options, statisticals &stats)	// impotiert eine SAM-Datei und speichert sie im Format candidates
{
	candidate tmp;	// tmp ist ein tempor�rer Speicher f�r ein 'candidate'.
	medianTree median;	// Es wird eine 'medianTree' ertzeugt f�r die Ermittelung des Medians.
	unsigned n=0;	//	Z�hler wird auf 0 gesetzt und z�hlt die Anzahl der Reads.
	
    // BamStream �ffnet einen Eingangsdatenstrom.
    seqan::BamStream SamInput(toCString(options.inputFile));
    if (!isGood(SamInput))
    {
		std::cerr << "-ERROR- Could not open" << std::endl;	// Ausgabe bei fehlerhafter Eingabe.
		return 1;
    }

	for(unsigned i=0; i<length(SamInput.header.sequenceInfos);i++)
	{
		seqan::appendValue(header,SamInput.header.sequenceInfos[i].i1);
	}
	// record enth�lt die Inhalte aus einer SAM-Datei.
	seqan::BamAlignmentRecord record;

	while (!atEnd(SamInput))
	{
		readRecord(record, SamInput);
		n++;	// Z�hler wird erh�ht.
		if(record.tLen>0)
		{
			median.add(record.tLen);	// Einf�gen eines neuen WErtes f�r die Medianberechnung.
			tmp.addvalues(record.rID,record.beginPos,record.tLen,candidate::undefine);	// tmp wird gef�llt.
			save.push_back(tmp);	// tmp wird abgespeichert.
		}
	}
	
	// Es wird der Median und die Standardabweichung gebildet mit R�ckfrage auf eine Eingabe durch die Kommandozeile.
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
	// Coverage wird berechnet.
	coverage(options,stats,length(record.seq),n,SamInput.header.sequenceInfos[0].i2);
	
    return 0;
}

int devide(candidates &input, candidates &result, modifier &options, statisticals &stats)	// sucht nach insertions und deletions
{
	// Es wird unterschieden ob ein 'candidate' eine Insertion oder eine Deletion ist.
	for (unsigned i=0;i<input.size();i++)	// Es wird �ber alle Elemente der Eingabe iteriert.
	{	// Fallunterscheidung:
		if(input[i].len<(stats.e_value-(options.x * stats.standardDegression)))	// 1. Fall : Deletion, da die Fragmentl�nge kleiner als Erwartungswert-x*Standardabweichung ist.
		{
			input[i].type=candidate::deletion;
			result.push_back(input[i]);
		}else
			if(input[i].len>(stats.e_value+(options.y * stats.standardDegression)))	// 2. Fall : Insertion, da die Fragmentl�nge gr��er als Erwartungswert-x*Standardabweichung ist.
			{
				input[i].type=candidate::insertion;
				result.push_back(input[i]);
			}
	}
	if(options.s){samSort(result);}	// Aufruf von samSort nach �berpr�fung der Eingabe der Kommandozeile.

	return 0;
}

int output(std::vector<std::vector<unsigned> > &input, seqan::String<seqan::CharString> &header, candidates &ref, char *out, statisticals &stats)
{
	std::fstream outFile(out, std::ios::out);	// Erzeugen eines Strems f�r den Output
	if (!outFile)
	{
		std::cerr << "-ERROR- Could not open output file" << std::endl;	// Ausgabe bei fehlerhafter funktion.
		return 1;
	}
	else
	{
		outFile << "ID\tAVEERAGE_POS\tAVEERAGE_LEN\tBegin\tEND\tTYPE\tCLUSTER_SIZE " << std::endl;	// Einf�gen einer Kopfzeile mit definitionen der Spalten.
		if(!input.empty())	// �berpr�fung ob Eingabae leer ist.
		{
			for (unsigned i=0;i<input.size();i++)	// Es wird �ber alle Elemente der Eingabe iteriert.
			{	// Initialisierung:
				unsigned targetSize=0;	// Speicher f�r die Durchschnitts L�nge.
				unsigned begin=ref[*input[i].begin()].pos;	// Speicher f�r den Anfang des Intervalls.
				unsigned end=0;	// Speicher f�r das Ende des Intervalls.
				unsigned n=input[i].size();	// Speicher f�r Clustergr��e.

				if(n>=stats.coverage)	// Grenzwert f�r die erzeugten Cluster vor der Ausgabe.
				{
					for(std::vector<unsigned>::iterator j=input[i].begin();j != input[i].end();j++)
					{
						targetSize += ref[*j].len;

						if(ref[*j].pos > end){end=ref[*j].pos;}
					}
				
					targetSize /= n;
					// Schreiben in die Ausgabedatei
					outFile << header[ref[*input[i].begin()].ref] << "\t" << (begin+end)/2<< "\t" << targetSize << "\t" << begin << "\t" << end << "\t" << (candidate::form)ref[*input[i].begin()].type << "\t" << n << std::endl;
				}
			}
		}
	}

	outFile.close();

	return 0;
}

void findCluster(candidates &input, std::vector<std::vector<unsigned> > &dest, candidate::form indel, statisticals &stats)
{
    std::vector<unsigned> currentCluster;
    unsigned currentClusterEndPos;
    bool emptyFlag = true;
	double coverage = stats.fold_sequence_coverage;

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
					if((double)overlapLen / std::min(currentClusterEndPos - input[currentCluster.back()].pos + 1, input[i].pos + input[i].len - 1) >  coverage){
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
