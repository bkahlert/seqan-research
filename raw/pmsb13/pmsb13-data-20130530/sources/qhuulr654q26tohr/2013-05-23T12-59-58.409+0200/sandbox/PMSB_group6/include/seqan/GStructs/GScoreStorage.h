//Autor:Jakob

#ifndef GINGER_STRUCTS_GSCORESTORAGE_H
#define GINGER_STRUCTS_GSCORESTORAGE_H

namespace seqan {
// score type
template<typename TScore>
/**
 * \brief Struct to contain the momentarily highest score and the related maxId, to determine the best match.
 */
struct GScoreStorage {
    TScore maxScore;/**<best score*/
    CharString maxId;/**<id of the sequence related to the best score*/
	bool scoreInitialized;/**<determines if maxScore is empty or not, e.g. if a comparison is needed*/
    GScoreStorage() {
		scoreInitialized=false;
    }
};
}  // namespace seqan

#endif  // GINGER_STRUCTS_GSCORESTORAGE_H 
