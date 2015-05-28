#include <iostream>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/align.h>
#include <seqan/seq_io.h>

using namespace seqan;
using namespace std;

struct MasterSequence{
    CharString id;
    Dna5String seq;
    // TODO store constructor parameter directly via :
    MasterSequence(CharString i, Dna5String s) :id(i),seq(s) {}
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


int main(int argc, char const ** argv)
{
    
    // Extract option values and print them.
	String<char> inPath;
	String<char> outPath;
	int threshold;
	if (readCommandlineParameters(inPath, outPath, threshold, argc, argv))
	    return 1;
	
	std::cout << "inPath  \t" << inPath << '\n'
	<< "outPath   \t" << outPath << '\n'
	<< "threshold\t" << threshold << '\n';
	
	
	
	seqan::CharString id;
	seqan::Dna5String seq;
	String<MasterSequence> masterSequences;
	
	seqan::SequenceStream inStream(toCString(inPath));
	if (!isGood(inStream))
	{
	    std::cerr << "ERROR: Could not open the input file.\n";
	    return 1;
	}
	seqan::SequenceStream outStream(toCString(outPath), seqan::SequenceStream::WRITE);
	if (!isGood(outStream))
	{
	    std::cerr << "ERROR: Could not open the output file.\n";
	    return 1;
	}
	
	
	Score<int> scoring(2, -1, -2, 0);
	LocalAlignmentEnumerator<Score<int>, Unbanded> enumerator(scoring, 5);
	Align< String<Dna> > align;
	resize(rows(align), 2);
	typedef Iterator<String<MasterSequence>, Rooted>::Type TMasterIt;
	TMasterIt masterIt;
	
	while (readRecord(id, seq, inStream)) {
	    assignSource(row(align, 0), seq);
	    bool createNewCluster = true;
	    masterIt = begin(masterSequences);
	    for (;!atEnd(masterIt);goNext(masterIt)){
		assignSource(row(align, 1), (*masterIt).seq);
		if (threshold <= localAlignment(align, scoring)){
		    createNewCluster=false;
		    append(id," $" + toCString((*masterIt).id) + "$");
		    writeRecord(id, seq, outStream);
		    break;
		}
	    }
	    if (createNewCluster){
		appendValue(masterSequences,MasterSequence(id, seq);
		append(id," $" + id + "$");
		writeRecord(id, seq, outStream);
	    }
	}	
	return 0;
}