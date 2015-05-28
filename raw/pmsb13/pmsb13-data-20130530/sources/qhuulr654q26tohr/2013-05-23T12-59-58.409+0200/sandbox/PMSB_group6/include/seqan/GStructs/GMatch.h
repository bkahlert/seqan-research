//Autor:Jakob

#ifndef GINGER_STRUCTS_GMATCH_H
#define GINGER_STRUCTS_GMATCH_H

namespace seqan {
// score type
template<typename TScore, typename TSequence>
/**
 * \brief Struct in which a computed match is saved
 */
struct GMatch {
	CharString id;/**<id of the target which is matched with the query*/
	CharString queryId;/**<id of the query*/
	TScore score;/**<score of the match*/
	TSequence querySequence;/**<sequence of the query*/
	TSequence targetSequence;/**<sequence of the target which is matched with the query*/
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
