#include <iostream>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/align.h>
#include <seqan/seq_io.h>

using namespace seqan;
using namespace std;

struct masterSequence{
    CharString id;
    Dna5String seq;
    // TODO store constructor parameter directly via :
    masterSequence(CharString i, Dna5String s) {
	id=i;
	s=seq;
    }
};


template<typename TPath, typename TThres>
int readCommandlineParameters(TPath & inPath, TPath & outPath, TThres & threshold, int & argc, char const ** & argv){
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("cluster");
    
    addArgument(parser, seqan::ArgParseArgument(
	seqan::ArgParseArgument::STRING, "Path to input fasta"));
    addArgument(parser, seqan::ArgParseArgument(
	seqan::ArgParseArgument::STRING, "Path to output fasta"));
    addOption(parser, seqan::ArgParseOption(
	"t", "threshold", "threshold for clustering local alignment score",
	seqan::ArgParseArgument::INTEGER, "INT"));
    
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    
    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
	return res == seqan::ArgumentParser::PARSE_ERROR;
    
    getArgumentValue(inPath, parser, 0);
    getArgumentValue(outPath, parser, 1);
    getOptionValue(threshold, parser, "threshold");
    return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char const ** argv)
{
    
    // Extract option values and print them.
	String<char> inPath;
	String<char> outPath;
	int threshold;
	if (readCommandlineParameters(inPath, outPath, threshold, argc, argv))
	    return 1;
	const char * constInPath = toCString(inPath);
	const char * constOutPath = toCString(outPath);
	
	std::cout << "inPath  \t" << inPath << '\n'
	<< "outPath   \t" << outPath << '\n'
	<< "threshold\t" << threshold << '\n';
	
	
	
	seqan::CharString id;
	seqan::Dna5String seq;
	StringSet<CharString> masterIds;
	
	seqan::SequenceStream inStream(constInPath);
	if (!isGood(inStream))
	{
	    std::cerr << "ERROR: Could not open the input file.\n";
	    return 1;
	}
	seqan::SequenceStream outStream(constOutPath, seqan::SequenceStream::WRITE);
	if (!isGood(outStream))
	{
	    std::cerr << "ERROR: Could not open the output file.\n";
	    return 1;
	}
	/*
	while (readRecord(id, seq, seqStream))
	std::cout << id << '\t' << seq << '\n';
	
	Align< String<Dna> > ali2;
	resize(rows(ali2), 2);
	assignSource(row(ali2, 0), "ataagcgtctcg");
	assignSource(row(ali2, 1), "tcatagagttgc");
	
	
	Score<int> scoring(2, -1, -2, 0);
	LocalAlignmentEnumerator<Score<int>, Unbanded> enumerator(scoring, 5);
	while (nextLocalAlignment(ali2, enumerator))
	{
	    std::cout << "Score = " << getScore(enumerator) << std::endl;
	    std::cout << ali2;
	    std::cout << "Aligns Seq1[" << clippedBeginPosition(row(ali2, 0)) << ":" << (clippedEndPosition(row(ali2, 0))-1) << "]";
	    std::cout << " and Seq2[" << clippedBeginPosition(row(ali2, 1)) << ":" <<  (clippedEndPosition(row(ali2, 1))-1) << "]" << std::endl << std::endl;
	}*/
	
	return 0;
}