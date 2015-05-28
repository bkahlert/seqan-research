//Autor:Jakob

#ifndef GINGER_STRUCTS_GMATCH_H
#define GINGER_STRUCTS_GMATCH_H

namespace seqan {
// score type
template<typename TScore, typename TSequence>
struct GMatch {
	CharString id;
	CharString queryId;
	TScore score;
	TSequence querySequence;
	TSequence targetSequence;
	GMatch(CharString i, TScore s) :
	id(i), score(s) {
	}
	GMatch() {
	}
	TScore getScore(){
		return score;
	}
};
}  // namespace seqan

#endif  // GINGER_STRUCTS_GMATCH_H 
