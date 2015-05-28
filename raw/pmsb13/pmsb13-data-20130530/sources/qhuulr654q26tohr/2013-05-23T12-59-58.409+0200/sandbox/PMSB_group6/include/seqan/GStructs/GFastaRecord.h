//Autor:Jakob

#ifndef GINGER_STRUCTS_GFASTARECORD_H
#define GINGER_STRUCTS_GFASTARECORD_H

namespace seqan {
	template<typename TSequence>
/**
* a struct that contains the id and the sequence of a record.
*/
	struct GFastaRecord {
		CharString id;/**<the id*/
		TSequence seq;/**<the sequence*/
		GFastaRecord(CharString i, TSequence s) :
		id(i), seq(s) {
		}
		GFastaRecord() {
		}
	};
}  // namespace seqan

#endif  // GINGER_STRUCTS_GFASTARECORD_H 
