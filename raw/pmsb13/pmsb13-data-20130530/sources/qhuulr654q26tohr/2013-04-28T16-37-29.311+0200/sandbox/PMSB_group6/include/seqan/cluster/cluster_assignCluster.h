#ifndef SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_CLUSTER_CLUSTER_ASSIGNCLUSTER_H_
#define SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_CLUSTER_CLUSTER_ASSIGNCLUSTER_H_

namespace seqan {
	template<typename TStream, typename TRecord>
	int assignToCluster(TStream &outStream, TRecord &newRecord,
						TRecord &clusterRecord) {
		CharString newId = newRecord.id;
		append(newId, " $");
		append(newId, clusterRecord.id);
		append(newId, "$");
		if (writeRecord(outStream, newId, newRecord.seq) != 0) {
			std::cerr << "ERROR: Could not write to file!\n";
			return 1;
		}
		return 0;
						}
}  // namespace seqan

#endif  // SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_CLUSTER_CLUSTER_ASSIGNCLUSTER_H_
