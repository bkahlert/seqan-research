#ifndef SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_CLUSTER_ALIGN_H_
#define SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_CLUSTER_ALIGN_H_

#include <seqan/align.h>

namespace seqan {
	template<typename TScore, typename TRecord, typename TScoringScheme>
	int myAlignFunc(TScore &score, TRecord &newRecord, TRecord &clusterRecord,
					TScoringScheme &scoringScheme) {
		Align<String<Dna5> > align;
		resize(rows(align), 2);
		assignSource(row(align, 0), newRecord.seq);
		assignSource(row(align, 1), clusterRecord.seq);
		score = localAlignment(align, scoringScheme);
		return 0;
					}
}  // namespace seqan

#endif  // SANDBOX_PMSB_GROUP6_INCLUDE_SEQAN_CLUSTER_ALIGN_H_
