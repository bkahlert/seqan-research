#ifndef GINGER_STRUCTS_GFASTARECORD_H
#define GINGER_STRUCTS_GFASTARECORD_H

namespace seqan {
	template<typename TSequence>
	struct GFastaRecord {
		CharString id;
		TSequence seq;
		GFastaRecord(CharString i, TSequence s) :
		id(i), seq(s) {
		}
		GFastaRecord() {
		}
	};
}  // namespace seqan

#endif  // GINGER_STRUCTS_GFASTARECORD_H 
