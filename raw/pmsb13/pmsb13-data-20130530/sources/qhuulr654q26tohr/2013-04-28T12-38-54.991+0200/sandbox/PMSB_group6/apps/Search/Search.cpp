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
#include <seqan/file.h>
#include <fstream>

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

template<typename TScore>
struct ScoredSequence {
      CharString id;
	  TScore score;
	  ScoredSequence(CharString i, TScore s) :
		id(i), score(s) {
	}
	ScoredSequence() {
	}
};


template<typename TPath>
int readCommandlineParameters(TPath & inCluster, TPath & inMaster, TPath & inQueries, TPath & outFile, int & argc, char const ** & argv) {
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
	
	

	// Parse command line.
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;

	getArgumentValue(inCluster, parser, 0);

	getArgumentValue(inMaster, parser, 1);

	getArgumentValue(inQueries, parser, 2);

	getArgumentValue(outFile, parser, 3);
	return 0;
}



int checkStreams(SequenceStream &myCluster, SequenceStream &myMaster, SequenceStream &myQueries, fstream &myOutput) {
	if (!isGood(myCluster)) {
		std::cerr << "ERROR: Could not open the Cluster file.\n";
		return 1;
	}

	if (!isGood(myMaster)) {
		std::cerr << "ERROR: Could not open the MasterSequences file.\n";
		return 1;
	}
	if (!isGood(myQueries)) {
		std::cerr << "ERROR: Could not open the Queries file.\n";
		return 1;
	}
	if (!myOutput.good()) {
		std::cerr << "ERROR: Could not open the output file.\n";
		return 1;   
	}
	return 0;
}

CharString getMasterId(String<char> & RefSeq){
  CharString MasterId;
  int dollarsFound=0; 
  for (int i=0;i<length(RefSeq);++i){
    if (RefSeq[i]=='$')
		++dollarsFound;
	if (dollarsFound==1)
      appendValue(MasterId,RefSeq[i]);
  }
  return MasterId;
}


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

	String<char> alpha;
	String<char> beta;
	String<char> inQueries;
	String<char> outFile;
	String<char> inMaster;
	String<char> inCluster;
	if (readCommandlineParameters(inCluster, inMaster, inQueries, outFile, argc, argv))
		return 1;
	SequenceStream myCluster(toCString(inCluster));
	SequenceStream myMaster(toCString(inMaster));
	SequenceStream myQueries(toCString(inQueries));
	std::fstream myOutput(toCString(outFile), std::ios::binary | std::ios::out);
	if (checkStreams(myCluster, myMaster, myQueries, myOutput))
		return 1;

	String<TRecord> masterSequences;
	TScoringScheme scoringScheme(1, -2, -1);
	TMasterIt masterIt;
	TRecord queryRecord;
	TRecord referenceRecord;
	TRecord clusterRecord;
	
	String<char> maxId;
	int score;
	ScoredSequence PartResults;
	String<ScoredSequence> Results;
	
	while (!readRecord(queryRecord.id, queryRecord.seq, myQueries)) {
		open(myMaster, toCString(inMaster));
		int maxScore=-1000;//nur provisorisch... muss überarbeitet werden
	      
	      while(!readRecord(referenceRecord.id, referenceRecord.seq, myMaster)){
			myAlignFunc(score, queryRecord, referenceRecord, scoringScheme);
			if (maxScore <= score) {
				maxScore=score;
				maxId=referenceRecord.id;
			}
	      }
	      while(!readRecord(clusterRecord.id, clusterRecord.seq, myCluster)){
		alpha=getMasterId(clusterRecord.id);
		beta=getMasterId(maxId);
		if (alpha==beta){
		  myAlignFunc(score, queryRecord, clusterRecord, scoringScheme);
		  stringstream ss << score;
		  PartResults.score = score.str();
		  PartResults.id=clusterRecord.id;
		  append(Results, PartResults);
		  }
		 }
	      /*sort(Results[].score); //muss noch gemacht werden */     
	     
	     
	      
	      seqan::CharString buffer;
	      append(buffer, "\n Query ID =");
	      append(buffer, queryRecord.id);
	     for(int i=0;i<length(Results);i++){
	      append(buffer, "\n \t MasterSequence ID =");
	      append(buffer, Results[i].id);
	      append(buffer, "\t AlignmentScore =");
	      append(buffer, Results[i].score);
	     }
	      
	      if(seqan::streamWriteBlock(myOutput, &buffer[0], length(buffer)) != length(buffer)){
		std::cerr << "ERROR: Could not print output.\n";
		return 1;}
	   }
	
	return 0;
}