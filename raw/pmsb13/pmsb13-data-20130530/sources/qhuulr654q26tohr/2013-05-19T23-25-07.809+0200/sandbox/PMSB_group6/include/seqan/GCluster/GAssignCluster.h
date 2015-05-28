#ifndef GINGER_CLUSTER_ASSIGNCLUSTER_H_
#define GINGER_CLUSTER_ASSIGNCLUSTER_H_

namespace seqan {
	template<typename TStream, typename TRecord>
	int assignToCluster(TStream &outStream, TRecord &newRecord,
						TRecord &clusterRecord) {
		CharString newId = newRecord.id;
		append(newId, "$");
		append(newId, clusterRecord.id);
		append(newId, "$");
		if (writeRecord(outStream, newId, newRecord.seq) != 0) {
			std::cerr << "ERROR: Could not write to file!\n";
			return 1;
		}
		return 0;
						}
}  // namespace seqan

#endif  // GINGER_CLUSTER_ASSIGNCLUSTER_H_
