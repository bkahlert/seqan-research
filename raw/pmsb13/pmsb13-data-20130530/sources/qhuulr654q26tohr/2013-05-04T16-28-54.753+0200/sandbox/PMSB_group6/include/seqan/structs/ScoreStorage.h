#ifndef SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_SCORESTORAGE_H
#define SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_SCORESTORAGE_H

namespace seqan {
	// score type
	template<typename TScore>
	struct ScoreStorage {
		TScore maxScore;
		CharString maxId;
		ScoreStorage() {
		}
	};
}  // namespace seqan

#endif  // SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_SCORESTORAGE_H 
