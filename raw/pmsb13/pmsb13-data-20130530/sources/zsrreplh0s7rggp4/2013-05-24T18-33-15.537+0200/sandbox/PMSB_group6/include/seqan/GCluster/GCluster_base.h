//Autor:Jakob

#ifndef GINGER_CLUSTER_CLUSTER_BASE_H_
#define GINGER_CLUSTER_CLUSTER_BASE_H_

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <seqan/GStructs/GFastaRecord.h>
#include <seqan/GStructs/PerformanceSample.h>
#include <seqan/GStructs/MemorySample.h>
#include <seqan/GSearch/GScore.h>
#include <seqan/GSearch/GScorePeptideQGram.h>
#include <seqan/GCluster/GAssignCluster.h>
#include <seqan/GCluster/GCheckClusterAssign.h>


using namespace seqan;
using namespace std;
/**
 * \brief Checks if the given streams are good.
 */
int checkStreams(SequenceStream &inStream/**<[in]the input stream*/, SequenceStream &outStream/**<[out]the clustered database output*/, SequenceStream &outMasterStream/**<[out]the mastersequences output*/) {
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
/**
 * \brief Searches for a sequence similar to a given query in a clustered database.
 * 
 * This is the base file of the GCluster algorithm. It includes all relevant functions and executes them.
 * Therefore the algorithm takes the first entry of a database and opens a new cluster with it. Every following entry is compared to all existing clusters to determine if it is added to that cluster or if it creates a new one.
 * If a sequence is added to a cluster, it gets the masterId of that cluster in addition to its own Id. If it opens a new cluster it gets is own Id as a master Id and therefore is the master of a new cluster.   
 */
int GCluster(String<PerformanceSample> & performance/**<[out]Struct in which the performance data is written*/, const char * inPath/**<[in]path for the data to be clustered*/, const char * outPath/**<[out]outpath for the clustered database*/, const char * outMasterPath/**<[out]outpath for the masterfile*/,
			 TSequence & sequenceFormatExample/**<[in]dummy object to define TSequence without defining other templates*/, TScore threshold/**<[in]threshold which is used to determine if a sequence creates a new cluster*/, double lengthThreshold/**<[in]threshold to determine how big the relative difference in length of two sequences may be so that thay are still compared*/, string & scoreMode/**<[in]the scoring mode which is used*/) {
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
	unsigned recordCounter=0;
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
			if (IsSameType<TSequence, Peptide>::VALUE){
				GScorePeptideQGram(score, record, *masterIt, scoringScheme, scoreMode);}
			else
				GScore(score, record, *masterIt, scoringScheme, scoreMode);
            perfGScoreAvg.end();
            appendValue(perfGScoreAvgContainer, perfGScoreAvg);
			if (GCheckClusterAssign(score, threshold, record, *masterIt, lengthThreshold)){
                createNewCluster = false;
				GAssignCluster(outStream, record, *masterIt);
                break;
            }
        }
        if (createNewCluster) {
            appendValue(masterSequences, TRecord(record.id, record.seq));
			GAssignCluster(outStream, record, record);
        }
        ++recordCounter;
        perfEntryAvg.end();
        appendValue(perfEntryAvgContainer, perfEntryAvg);
    }
    for (masterIt=begin(masterSequences); !atEnd(masterIt); goNext(masterIt)) {
        writeRecord(outMasterStream, (*masterIt).id, (*masterIt).seq);
	}
	if (length(masterSequences)>0)
		cout << "\tCreated "<< length(masterSequences) << " cluster with an average of " << ((double)recordCounter/length(masterSequences)) << " sequences each."<< endl;
	else
		cout << "\tCreated 0 cluster." << endl;
	
    //MemorySample ms("End of cluster function");
    //ms.printToStdout();

    PerformanceSample perfEntryAvg("Handle one database entry in Cluster average");
    //perfEntryAvg.takeAverageValues(perfEntryAvgContainer);
    //appendValue(performance, perfEntryAvg);

    PerformanceSample perfGScoreAvg("Call GScore in Cluster average");
    //perfGScoreAvg.takeAverageValues(perfGScoreAvgContainer);
    //appendValue(performance, perfGScoreAvg);

    perfCluster.end();
    appendValue(performance, perfCluster);
    return 0;
}

#endif  // GINGER_CLUSTER_CLUSTER_BASE_H_