#ifndef GINGER_STRUCTS_GMATCH_H
#define GINGER_STRUCTS_GMATCH_H

namespace seqan {
// score type
template<typename TScore>
struct GMatch {
	CharString id;
	TScore score;
	GMatch(CharString i, TScore s) :
	id(i), score(s) {
	}
	GMatch() {
	}
};
}  // namespace seqan

#endif  // GINGER_STRUCTS_GMATCH_H 
