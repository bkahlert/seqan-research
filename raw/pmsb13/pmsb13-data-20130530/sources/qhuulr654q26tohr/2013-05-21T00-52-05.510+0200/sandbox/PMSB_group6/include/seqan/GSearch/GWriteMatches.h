//Autor:Hannes

#ifndef GINGER_GSEARCH_GWRITEMATCHES_H_
#define GINGER_GSEARCH_GWRITEMATCHES_H_

#include <seqan/sequence.h>
#include <fstream>
#include <seqan/GStructs/GMatch.h>
#include <seqan/GStructs/GFastaRecord.h>
#include <seqan/getIdWithoutMaster.h>

using namespace seqan;
using namespace std;

template<typename TScore, typename TSequence>
int GWriteMatches(std::fstream & outputStream, String<GMatch<TScore,TSequence> > & matches, GFastaRecord<TSequence> & queryRecord)
{
	stringstream ss;
	seqan::CharString buffer;
	
	
	if(outputStream.tellp()==0){
	ss.str("");
	ss.clear();
	ss << "QueryRecordID\t";
	ss << "MatchId\t";
	ss << "Q-Gram-Score\t";
	ss << "Global-Alignment-Score\n";
	append(buffer, ss.str());	  
	}
	
	for(int i=0; i<length(matches); i++) {
		Align <TSequence> align;
		resize(rows(align), 2);
		assignSource(row(align, 0), queryRecord.seq);
		assignSource(row(align, 1), matches[i].targetSequence);
		int tempScore = globalAlignment(align,Score<int,Simple>(1,-1,-1),AlignConfig<true, true, true, true>());
	
		ss.str("");
		ss.clear();
		ss << queryRecord.id;
		ss << "\t";
		ss << getIdWithoutMaster(matches[i].id);
		ss << "\t";
		ss << matches[i].score;
		ss << "\t";
		ss << tempScore;
		ss << "\n";
		//ss << align;
		//ss << "\n";
		append(buffer, ss.str());
		
		}
		outputStream.write(toCString(buffer), length(buffer));
		/*
	if() {
		std::cerr << "ERROR: Could not print output.\n";
		return 1;
	}*/
	clear(matches);
	return 0;
}

#endif  // GINGER_GSEARCH_GWRITEMATCHES_H_