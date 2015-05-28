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
	
	append(buffer, "Query ID =");
	append(buffer, queryRecord.id);
	
	for(int i=0; i<length(matches); i++) {
	
		Align < String<Dna5> > align;
		resize(rows(align), 2);
		assignSource(row(align, 0), queryRecord.seq);
		assignSource(row(align, 1), matches[i].targetSequence);
		int tempScore = globalAlignment(align,Score<int,Simple>(1,-1,0),Hirschberg());
	
		append(buffer, queryRecord.id);
		append(buffer, "\t");
		append(buffer, getIdWithoutMaster(matches[i].id));
		append(buffer, "\t");
		ss.str("");
		ss.clear();
		ss << matches[i].score;
		append(buffer, ss.str());
		append(buffer, "\t");
		append(buffer, length(align));
		append(buffer, "\t");
		ss.str("");
		ss.clear();
		ss << tempScore;
		append(buffer, ss.str());
		}
	append(buffer, "\n");
	if(seqan::streamWriteBlock(outputStream, &buffer[0], length(buffer)) != length(buffer)) {
		std::cerr << "ERROR: Could not print output.\n";
		return 1;
	}
	clear(matches);
	return 0;
}

#endif  // GINGER_GSEARCH_GWRITEMATCHES_H_