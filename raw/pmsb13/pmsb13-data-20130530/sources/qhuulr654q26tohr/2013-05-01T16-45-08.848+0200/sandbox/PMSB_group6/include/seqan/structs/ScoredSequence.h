#ifndef SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_SCOREDSEQUENCE_H
#define SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_SCOREDSEQUENCE_H

namespace seqan {
	struct ScoredSequence {
		CharString id;
		String<int> score;
		ScoredSequence(CharString i, String<int> s) :
		id(i), score(s) {
		}
		ScoredSequence() {
		}
	};
}  // namespace seqan

#endif  // SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_SCOREDSEQUENCE_H 
