//autor:Jakob

#ifndef GINGER_GCLUSTER_GCHECKCLUSTERASSIGN_H_
#define GINGER_GCLUSTER_GCHECKCLUSTERASSIGN_H_

#include <seqan/GStructs/GFastaRecord.h>

namespace seqan {
	template<typename TScore, typename TSequence>
	/** \brief checks whether the new record should be assigned to the cluster (represented by clusterRecord).
	 * The function checks whether the new record should be assigned to the cluster (represented by clusterRecord).
	 * If the given score is greater or equal the threshold and if the length of the new record's sequence does not exceed (<= and >=) the length of the cluster master (+- relative lengthThreshold), the function returns true.
	 */
	bool GCheckClusterAssign(TScore &score/**<[in] score of comparison of newRecord to clusterMaster (clusterRecord)*/, TScore &threshold/**<[in] threshold to compare against*/, GFastaRecord<TSequence> &newRecord/**<[in] record which is either assigned to a cluster or forms a new cluster*/, GFastaRecord<TSequence> &clusterRecord/**<[in] cluster master of the cluster to which the newRecord would be assigned*/, double &lengthThreshold/**<[in] relative acceptable sequence length difference of cluster master to newRecord*/) {
		return (score>=threshold&&length(newRecord.seq)<=(1.0+lengthThreshold)*length(clusterRecord.seq)&&length(newRecord.seq)>=(1.0-lengthThreshold)*length(clusterRecord.seq));
	}
}  // namespace seqan

#endif  // GINGER_GCLUSTER_GCHECKCLUSTERASSIGN_H_
