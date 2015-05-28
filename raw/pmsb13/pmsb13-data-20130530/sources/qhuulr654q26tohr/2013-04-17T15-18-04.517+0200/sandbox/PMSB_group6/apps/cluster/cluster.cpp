#include <iostream>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/align.h>
#include <seqan/seq_io.h>

using namespace seqan;
using namespace std;

/*struct MasterSequence {
	CharString id;
	Dna5String seq;
	MasterSequence(CharString i, Dna5String s) :
		id(i), seq(s) {
	}
};*/

template<typename TSequence>
struct FastaRecord {
	CharString id;
	TSequence seq;
	FastaRecord(CharString i, TSequence s) :
		id(i), seq(s) {
	}
	FastaRecord() {
	}
};

template<typename TPath, typename TThres>
int readCommandlineParameters(TPath & inPath, TPath & outPath,
		TThres & threshold, int & argc, char const ** & argv) {
	// Setup ArgumentParser.
	seqan::ArgumentParser parser("cluster");

	addArgument(
			parser,
			seqan::ArgParseArgument(seqan::ArgParseArgument::STRING,
					"Path to input fasta"));
	addArgument(
			parser,
			seqan::ArgParseArgument(seqan::ArgParseArgument::STRING,
					"Path to output fasta"));
	addOption(
			parser,
			seqan::ArgParseOption("t", "threshold",
					"threshold for clustering local alignment score",
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

int checkStreams(SequenceStream &inStream, SequenceStream &outStream){
    if (!isGood(inStream)) {
	std::cerr << "ERROR: Could not open the input file.\n";
	return 1;
    }
    if (!isGood(outStream)) {
	std::cerr << "ERROR: Could not open the output file.\n";
	return 1;
    }
    return 0;
}

template<typename TScore, typename TRecord, typename TScoringScheme>
int myAlignFunc(TScore &score, TRecord &newRecord, TRecord &clusterRecord, TScoringScheme &scoringScheme){
	Align < String<Dna5> > align;
	resize(rows(align), 2);
	assignSource(row(align, 0), newRecord.seq);
	assignSource(row(align, 1), clusterRecord.seq);
	score = localAlignment(align, scoringScheme);
	return 0;
}

template<typename TStream, typename TRecord>
assignToCluster(TStream outStream, TRecord &newRecord, TRecord &clusterRecord){
    CharString newId = newRecord.id;
    append(newId, " $");
    append(newId, clusterRecord.id);
    append(newId, "$");
    if (writeRecord(outStream, newId, newRecord.seq) != 0) {
	std::cerr << "ERROR: Could not write to file!\n";
	return 1;
    }
    return 0;
}


int main(int argc, char const ** argv) {
    typedef int TThreshold;
    typedef Score<int> TScoringScheme;
    typedef FastaRecord<Dna5String> TRecord;
    typedef Iterator<String<TRecordString> > , Rooted>::Type TMasterIt;
    typedef int TScore;
    
	String<char> inPath;
	String<char> outPath;
	TThreshold threshold;
	if (readCommandlineParameters(inPath, outPath, threshold, argc, argv))
	    return 1;
	SequenceStream inStream(toCString(inPath));
	SequenceStream outStream(toCString(outPath), SequenceStream::WRITE);
	if (checkStreams(inStream, outStream))
	    return 1;
	
	String<TRecord> masterSequences;
	TScoringScheme scoringScheme(2, -1, -2, 0);
	TMasterIt masterIt;
	TRecord record;
	
	while (!readRecord(record.id, record.seq, inStream)) {
		bool createNewCluster = true;
		masterIt = begin(masterSequences);
		for (; !atEnd(masterIt); goNext(masterIt)) {
			TScore score;
			myAlignFunc(score, record, *masterIt, scoringScheme);
			if (threshold <= score) {
				createNewCluster = false;
				assignToCluster(outStream, record, *masterIt);
				break;
			}
		}
		if (createNewCluster) {
		    appendValue(masterSequences, MasterSequence(record.id, record.seq));
		    assignToCluster(outStream, record, record);
		}
	}
	return 0;
}