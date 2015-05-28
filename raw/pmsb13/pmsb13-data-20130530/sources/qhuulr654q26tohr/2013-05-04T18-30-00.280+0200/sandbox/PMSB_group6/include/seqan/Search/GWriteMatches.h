#ifndef GINGER_SEARCH_GWRITEMATCHES_H_
#define GINGER_SEARCH_GWRITEMATCHES_H_

#include <seqan/sequence.h>
#include <fstream>
#include <seqan/structs/GMatch.h>
#include <seqan/structs/GFastaRecord.h>

using namespace seqan;
using namespace std;

template<typename TScore, typename TSequence>
int GWriteMatches(std::fstream & outputStream, String<GMatch<TScore> > & matches, GFastaRecord<TSequence> & queryRecord)
{
	stringstream ss;
	seqan::CharString buffer;
	
	append(buffer, "\n Query ID =");
	append(buffer, queryRecord.id);
	
	for(int i=0; i<length(matches); i++) {
		append(buffer, "\n \t MasterSequence ID =");
		append(buffer, matches[i].id);
		append(buffer, "\t AlignmentScore =");
		ss.str("");
		ss.clear();
		ss << matches[i].score;
		append(buffer, ss.str());
	}
	if(seqan::streamWriteBlock(outputStream, &buffer[0], length(buffer)) != length(buffer)) {
		std::cerr << "ERROR: Could not print output.\n";
		return 1;
	}
	return 0;
}

#endif  // GINGER_SEARCH_GWRITEMATCHES_H_