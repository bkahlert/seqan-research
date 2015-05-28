// ==========================================================================
//                                   Search
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Hannes Blumenthal <hannes.blumenthal@gmx.net>
// ==========================================================================
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/align.h>
#include <seqan/seq_io.h>

using namespace seqan;
using namespace std;

/*
2.:Globales Alignment der gesuchten sequenzen in der Datenbank
3.: Ausgeben der gesuchten sequenzen mit ähnlichster datenbanksequenz und ähnlichkeit
*/

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
int readCommandlineParameters(TPath & inCluster, TPath & inMaster, TPath & inQueries, TPath & outFile,
		TThres & threshold, int & argc, char const ** & argv) {
	// Setup ArgumentParser.
	seqan::ArgumentParser parser("search");

	addArgument(
			parser,
			seqan::ArgParseArgument(seqan::ArgParseArgument::STRING,
					"Path to input cluster FASTA"));
	addArgument(
			parser,
			seqan::ArgParseArgument(seqan::ArgParseArgument::STRING,
					"Path to input Master FASTA"));
	addArgument(
			parser,
			seqan::ArgParseArgument(seqan::ArgParseArgument::STRING,
					"Path to input queries FASTA"));
	addArgument(
			parser,
			seqan::ArgParseArgument(seqan::ArgParseArgument::STRING,
					"Path to output file FASTA"));
	
	addOption(
			parser,
			seqan::ArgParseOption("t", "threshold",
					"threshold for searching global alignment score",
					seqan::ArgParseArgument::INTEGER, "INT"));

	// Parse command line.
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;

	getArgumentValue(inQueries, parser, 0);

	getArgumentValue(inCluster, parser, 1);

	getArgumentValue(inMaster, parser, 2);

	getArgumentValue(outFile, parser, 3);
	getOptionValue(threshold, parser, "threshold");
	return 0;
}



int checkStreams(SequenceStream &myQueries, SequenceStream &myCluster, SequenceStream &myMaster, SequenceStream &myOutput) {
	if (!isGood(myQueries)) {
		std::cerr << "ERROR: Could not open the queries file.\n";
		return 1;
	}
	if (!isGood(myOutput)) {
		std::cerr << "ERROR: Could not open the output file.\n";
		return 1;
	}
	if (!isGood(myMaster)) {
		std::cerr << "ERROR: Could not open the MasterSequences file.\n";
		return 1;
	}
	if (!isGood(myCluster)) {
		std::cerr << "ERROR: Could not open the Cluster file.\n";
		return 1;
	}
	return 0;
}


/*int main(int argc, char const ** argv)
{
    if (argc < 4)
    {
        std::cerr << "USAGE: basic_seq_io_example FILENAME\n";
        return 1;
    }
    seqan::CharString id;
    seqan::Dna5String seq;
    seqan::SequenceStream cluster(argv[1]);
    seqan::SequenceStream clusterMaster(argv[2]);
    seqan::SequenceStream queries(argv[3]);
    if (!isGood(file))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }
    while (!atEnd(file))
    {
        if (readRecord(id, seq, seqStream) != 0)
        {
            std::cerr << "ERROR: Could not read from example.fa!\n";
            return 1;
        }
        std::cout << id << '\t' << seq << '\n';
    }
    return 0;
}*/



template<typename TScore, typename TRecord, typename TScoringScheme>
int myAlignFunc(TScore &score, TRecord &newRecord, TRecord &clusterRecord,
		TScoringScheme &scoringScheme) {
	Align < String<Dna5> > align;
	resize(rows(align), 2);
	assignSource(row(align, 0), newRecord.seq);
	assignSource(row(align, 1), clusterRecord.seq);
	score = globalAlignment(align, scoringScheme);
	return 0;
}

int main(int argc, char const ** argv) {
	typedef int TThreshold;
	typedef Score<int> TScoringScheme;
	typedef FastaRecord<Dna5String> TRecord;
	typedef Iterator<String<TRecord> , Rooted>::Type TMasterIt;
	typedef int TScore;

	String<char> inQueries;
	String<char> outFile;
	String<char> inMaster;
	String<char> inCluster;
	if (readCommandlineParameters(inCluster, inMaster, inQueries, outFile, argc, argv))
		return 1;
	SequenceStream myCluster(toCString(inCluster));
	SequenceStream myMaster(toCString(inMaster));
	SequenceStream myQueries(toCString(inQueries));
	SequenceStream myOutput(toCString(outFile), SequenceStream::WRITE);
	if (checkStreams(inCluster, inMaster, inQueries, outFile))
		return 1;

	String<TRecord> masterSequences;
	TScoringScheme scoringScheme(1, -2, -1);
	TMasterIt masterIt;
	TRecord queryRecord;
	TRecord referenceRecord;
	maxScore;
	maxId;

	while (readRecord(queryRecord.id, queryRecord.seq, myQueries)) {
		open(myMaster, inMaster);
	      while(readRecord(referenceRecord.id, referenceRecord.seq, myMaster)){
			myAlignFunc(score, queryRecord, referenceRecord, scoringScheme);
			if (maxScore <= score) {
				maxScore=score;
				maxId=referenceRecord.id;
			}
	  
	      seqan::CharString buffer;
	      append(buffer, "\n Query ID =");
	      append(buffer, queryRecord.id);
	      append(buffer, "\t MasterSequence ID =");
	      append(buffer, maxId);
	      append(buffer, "\t AlignmentScore =");
	      append(buffer, maxScore);
	      
	      seqan::streamWriteBlock(myOutput, &buffer[0], length(buffer)) != length(buffer);
	   }
	
	return 0;
}