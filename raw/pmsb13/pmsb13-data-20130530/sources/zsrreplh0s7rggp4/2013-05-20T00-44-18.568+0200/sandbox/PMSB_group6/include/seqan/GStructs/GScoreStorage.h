//Autor:Jakob

#ifndef GINGER_STRUCTS_GSCORESTORAGE_H
#define GINGER_STRUCTS_GSCORESTORAGE_H

namespace seqan {
// score type
template<typename TScore>
struct GScoreStorage {
    TScore maxScore;
    CharString maxId;
	bool scoreInitialized;
    GScoreStorage() {
		scoreInitialized=false;
    }
};
}  // namespace seqan

#endif  // GINGER_STRUCTS_GSCORESTORAGE_H 
