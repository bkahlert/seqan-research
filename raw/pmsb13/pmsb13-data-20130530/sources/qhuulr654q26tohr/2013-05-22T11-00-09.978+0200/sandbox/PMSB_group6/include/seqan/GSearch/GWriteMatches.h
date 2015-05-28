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
/** \brief writes the given matches into an fstream
 *  
 *  This function gets the reference to an outputStream(fstream), a string of GMatch and a queryRecord(GFastaRecord).
 *  It computes a global alignment without overlap penalaties for every found hit. 
 *  In the outputStream it writes the queryRecordId, the Id of the Matches, a saved Q-Gram Score and a global alignment score.
 */
int GWriteMatches(std::fstream & outputStream /**<[out] stream in which is written*/, String<GMatch<TScore,TSequence> > & matches/**<[in] String of Matches which will be written to the stream*/, GFastaRecord<TSequence> & queryRecord/**<[in] queryRecord from wich the queryRecordId is taken to compute the global alignment*/)
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
	clear(matches);
	return 0;
}

#endif  // GINGER_GSEARCH_GWRITEMATCHES_H_