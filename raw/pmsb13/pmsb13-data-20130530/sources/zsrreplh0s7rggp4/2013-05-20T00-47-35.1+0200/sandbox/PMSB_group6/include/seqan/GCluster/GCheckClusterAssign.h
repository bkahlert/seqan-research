#ifndef GINGER_GCLUSTER_GCHECKCLUSTERASSIGN_H_
#define GINGER_GCLUSTER_GCHECKCLUSTERASSIGN_H_

#include <seqan/GStructs/GFastaRecord.h>

namespace seqan {
	template<typename TScore, typename TSequence>
	bool GCheckClusterAssign(TScore &score, TScore &threshold, GFastaRecord<TSequence> &newRecord, GFastaRecord<TSequence> &clusterRecord, double &lengthThreshold) {
		return (score>=threshold&&length(newRecord)<(1.0+lengthThreshold)*length(clusterRecord)&&length(newRecord)>(1.0-lengthThreshold)*length(clusterRecord));
	}
}  // namespace seqan

#endif  // GINGER_GCLUSTER_GCHECKCLUSTERASSIGN_H_
