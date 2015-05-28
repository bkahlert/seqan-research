#ifndef GINGER_GCLUSTER_GCHECKCLUSTERASSIGN_H_
#define GINGER_GCLUSTER_GCHECKCLUSTERASSIGN_H_

#include <seqan/GStructs/GFastaRecord.h>

namespace seqan {
	template<typename TScore, typename TSequence>
	bool GCheckClusterAssign(TScore &score, TScore &threshold, GFastaRecord<TSequence> &newRecord, GFastaRecord<TSequence> &clusterRecord, double &lengthThreshold) {
		return (score>=threshold&&length(newRecord.seq)<(1.0+lengthThreshold)*length(clusterRecord.seq)&&length(newRecord.seq)>(1.0-lengthThreshold)*length(clusterRecord.seq));
	}
}  // namespace seqan

#endif  // GINGER_GCLUSTER_GCHECKCLUSTERASSIGN_H_
