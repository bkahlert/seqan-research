#ifndef SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_FASTARECORD_H
#define SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_FASTARECORD_H

namespace seqan {
	template<typename TSequence>
	struct FastaRecord {
		CharString id;
		TSequence seq;
		FastaRecord(CharString i, TSequence s) :
		id(i), seq(s) {
		}
		FastaRecord() {
		}
	};
}  // namespace seqan

#endif  // SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_STRUCTS_FASTARECORD_H 
