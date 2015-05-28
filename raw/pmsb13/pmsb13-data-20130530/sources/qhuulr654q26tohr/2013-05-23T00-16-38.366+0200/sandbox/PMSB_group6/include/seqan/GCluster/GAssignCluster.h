//Autor:Jakob

#ifndef GINGER_CLUSTER_ASSIGNCLUSTER_H_
#define GINGER_CLUSTER_ASSIGNCLUSTER_H_

namespace seqan {
	template<typename TStream, typename TRecord>
	/**
	 * \brief assigns a cluster ID to the relevant RecordID
	 * 
	 * this function gets a stream, in which the output is written, the record which should be assigned and 
	 * the cluster to which the record should be assigned, so that the output stream contains the record with a clusterID.
	 */
	int GAssignCluster(TStream &outStream/**<[out] stream, in which the output is written */, TRecord &newRecord/**<[in] record which will be assigned to a cluster */,
						TRecord &clusterRecord /**<[in] cluster to which the record will be assigned  */) {
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
