#ifndef SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_SCOREDSEQUENCE_H
#define SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_SCOREDSEQUENCE_H

namespace seqan {
	struct ScoredSequence {
		CharString id;
		CharString score;
		ScoredSequence(CharString i, CharString s) :
		id(i), score(s) {
		}
		ScoredSequence() {
		}
	};
}  // namespace seqan

#endif  // SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_SCOREDSEQUENCE_H 
