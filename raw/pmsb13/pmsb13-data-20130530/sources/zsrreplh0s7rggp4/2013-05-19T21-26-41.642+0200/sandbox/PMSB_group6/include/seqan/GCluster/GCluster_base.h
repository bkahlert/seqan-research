#ifndef GINGER_CLUSTER_CLUSTER_BASE_H_
#define GINGER_CLUSTER_CLUSTER_BASE_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <seqan/GStructs/GFastaRecord.h>
#include <seqan/GStructs/PerformanceSample.h>
#include <seqan/GStructs/MemorySample.h>
#include <seqan/GSearch/GScore.h>
#include <seqan/GCluster/GAssignCluster.h>
#include <seqan/GCluster/GCheckClusterAssign.h>


using namespace seqan;
using namespace std;

int checkStreams(SequenceStream &inStream, SequenceStream &outStream, SequenceStream &outMasterStream) {
    if (!isGood(inStream)) {
        std::cerr << "ERROR: Could not open the input file.\n";
        return 1;
    }
    if (!isGood(outStream)) {
        std::cerr << "ERROR: Could not open the output file.\n";
        return 1;
    }
    if (!isGood(outMasterStream)) {
        std::cerr << "ERROR: Could not open the output file.\n";
        return 1;
    }
    return 0;
}

template<typename TSequence, typename TScore>
int GCluster(String<PerformanceSample> & performance, const char * inPath, const char * outPath, const char * & outMasterPath,
			 TSequence & sequenceFormatExample, TScore threshold, double lengthThreshold, string & scoreMode) {
    PerformanceSample perfCluster("GCluster Total");
	typedef Score<int> TScoringScheme;
	typedef GFastaRecord<TSequence> TRecord;
    typedef typename Iterator<String<TRecord> , Rooted>::Type TMasterIt;
	
	SequenceStream inStream(inPath);
	SequenceStream outStream(outPath, SequenceStream::WRITE);
	SequenceStream outMasterStream(outMasterPath, SequenceStream::WRITE);
    if (checkStreams(inStream, outStream, outMasterStream))
        return 1;

    String<TRecord> masterSequences;
    TScoringScheme scoringScheme(1, -2, -1);
    TMasterIt masterIt;
    TRecord record;

    String<PerformanceSample> perfEntryAvgContainer;
    String<PerformanceSample> perfGScoreAvgContainer;
	
    while (!readRecord(record.id, record.seq, inStream)) {
        PerformanceSample perfEntryAvg("");
		
        bool createNewCluster = true;
        masterIt = begin(masterSequences);
        for (; !atEnd(masterIt); goNext(masterIt)) {
            TScore score;
            PerformanceSample perfGScoreAvg("");
            GScore(score, record, *masterIt, scoringScheme, scoreMode);
            perfGScoreAvg.end();
            appendValue(perfGScoreAvgContainer, perfGScoreAvg);
			if (GCheckClusterAssign(score, threshold, record, *masterIt, lengthThreshold)){
                createNewCluster = false;
                assignToCluster(outStream, record, *masterIt);
                break;
            }
        }
        if (createNewCluster) {
            appendValue(masterSequences, TRecord(record.id, record.seq));
            assignToCluster(outStream, record, record);
        }
        perfEntryAvg.end();
        appendValue(perfEntryAvgContainer, perfEntryAvg);
    }
    for (masterIt=begin(masterSequences); !atEnd(masterIt); goNext(masterIt)) {
        writeRecord(outMasterStream, (*masterIt).id, (*masterIt).seq);
    }

    //MemorySample ms("End of cluster function");
    //ms.printToStdout();

    PerformanceSample perfEntryAvg("Handle one database entry in Cluster average");
    perfEntryAvg.takeAverageValues(perfEntryAvgContainer);
    appendValue(performance, perfEntryAvg);

    PerformanceSample perfGScoreAvg("Call GScore in Cluster average");
    perfGScoreAvg.takeAverageValues(perfGScoreAvgContainer);
    appendValue(performance, perfGScoreAvg);

    perfCluster.end();
    appendValue(performance, perfCluster);
    return 0;
}

#endif  // GINGER_CLUSTER_CLUSTER_BASE_H_