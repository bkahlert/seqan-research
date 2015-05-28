/*
 *   Flexbar.h
 *
 *   Author: jtr
 */

#ifndef FLEXBAR_FLEXBAR_H_
#define FLEXBAR_FLEXBAR_H_

#include <string>
#include <iostream>

#include <tbb/pipeline.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_vector.h>

#include <seqan/basic.h>
#include <seqan/arg_parse.h>

#include "Enums.h"
#include "Options.h"
#include "AdapterLoader.h"
#include "SequencingRead.h"
#include "SequenceInputFilter.h"
#include "MultiplexedInputFilter.h"
#include "MultiplexedOutputFilter.h"
#include "MultiplexedAlignmentFilter.h"


void loadBarcodesAndAdapters(Options &o){
	
	using namespace std;
	using namespace flexbar;
	
	using seqan::CharString;
	
	if(o.barDetect != BOFF){
		
		tbb::task_scheduler_init init_serial(1);
		tbb::pipeline bpipeline;
		
		SequenceInputFilter<CharString, CharString > adapter_filter(o, o.barcodeFile, true, false);
		bpipeline.add_filter(adapter_filter);
		
		AdapterLoader<CharString, CharString> adapterLoader(o.format);
		bpipeline.add_filter(adapterLoader);
		
		bpipeline.run(1);
		o.barcodes = adapterLoader.getAdapters();
		
		adapterLoader.printAdapters("Barcode");
		
		if(o.barcodes.size() == 0){
			cerr << "No barcodes found in file!\n" << endl;
			exit(1);
		}
	}
	
	if(o.adapRm != AOFF){
		
		AdapterLoader<CharString, CharString> adapterLoader(o.format);
		
		if(o.useAdapterFile){
			tbb::task_scheduler_init init_serial(1);
			tbb::pipeline prepipeline;
			
			SequenceInputFilter<CharString, CharString> adapter_filter(o, o.adapterFile, true, false);
			prepipeline.add_filter(adapter_filter);
			
			prepipeline.add_filter(adapterLoader);
			prepipeline.run(1);
			
			o.adapters = adapterLoader.getAdapters();
			
			if(o.adapters.size() == 0){
				cerr << "No adapters found in file!\n" << endl;
				exit(1);
			}
		}
		else {
			SequencingRead<CharString, CharString> *myRead;
			myRead = new SequencingRead<CharString, CharString>(o.adapterSeq, "cmdline");
			
			TAdapter adap;
			adap.first = myRead;
			o.adapters.push_back(adap);
			
			adapterLoader.setAdapters(o.adapters);
		}
		
		adapterLoader.printAdapters("Adapter");
	}
}


void printComputationTime(const time_t t_start){
	using namespace std;
	
	time_t t_end;
	time(&t_end);
	
	int totalTime = int(difftime(t_end, t_start));
	
	int hours   = div(totalTime, 3600).quot;
	int rest    = div(totalTime, 3600).rem;
	int minutes = div(rest, 60).quot;
	int seconds = div(rest, 60).rem;
	
	cout << "Computation time:  ";
	if(hours > 0)                               cout << hours     << " h ";
	if(hours > 0 || minutes > 0)                cout << minutes   << " min ";
	if(hours > 0 || minutes > 0 || seconds > 0) cout << seconds   << " sec\n\n\n";
	else                                        cout              << "< 1 sec\n\n\n";
}


void startComputation(Options &o){
	
	using namespace std;
	using namespace flexbar;
	
	typedef seqan::CharString TString;
	typedef seqan::CharString TIDString;
	
	time_t t_start;
	time(&t_start);
	
	tbb::task_scheduler_init init_serial(o.nThreads);
	tbb::pipeline pipe;
	
	MultiplexedInputFilter<TString, TIDString > inputFilter(o);
	
	if(o.adapRm == AOFF && o.barDetect == BOFF) cout << "\n";
	cout << "\nProcessing reads ..." << flush;
	if(o.logLevel != NONE) cout << "\n\nGenerating " << o.log_level << " verbose output.\n" << endl;
	
	MultiplexedAlignmentFilter<TString, TIDString> alignFilter(o, &o.barcodes, &o.adapters);
	MultiplexedOutputFilter<TString, TIDString> outputFilter(o, &o.barcodes, &o.adapters);
	
	pipe.add_filter(inputFilter);
	pipe.add_filter(alignFilter);
	pipe.add_filter(outputFilter);
	
	pipe.run(o.nThreads);
	
	if(o.logLevel == TAB) cout << "\n";
	cout << "done.\n" << endl;
	
	printComputationTime(t_start);
	
	
	// barcode and adapter removal statistics
	
	if(o.writeLengthDist) outputFilter.writeLengthDist();
	
	if(o.adapRm != AOFF){
		outputFilter.printAdapterRemovalStats();
		alignFilter.printAdapterOverlapStats();
	}

	outputFilter.printFileSummary();


	const unsigned long nReads     = inputFilter.getNrProcessedReads();
	const unsigned long nGoodReads = outputFilter.getNrGoodReads();
	const unsigned long uncalled   = inputFilter.getNrUncalledReads();
	const unsigned long uPairs     = inputFilter.getNrUncalledPairedReads();

	cout << "Filtering statistics\n";
	cout << "====================\n";
	cout << "Processed reads                   " << nReads << endl;
	cout << "  skipped due to uncalled bases   ";

	if(o.runType == PAIRED || o.runType == PAIRED_BARCODED){
		cout << 2 * uPairs;
	
		if(uncalled > 0)
		cout << "   (" << uncalled << " uncalled in " << uPairs << " pairs)";
		cout << endl;
	}
	else cout << uncalled << endl;

	if(o.phred_preQual > 0)
	cout << "  trimmed due to low quality      " << inputFilter.getNrLowPhredReads() << endl;

	if(o.adapRm != AOFF)
	cout << "  short prior adapter removal     " << alignFilter.getNrPreShortReads() << endl;

	cout << "  finally skipped short reads     " << outputFilter.getNrShortReads() << endl;
	cout << "Discarded reads overall           " << (nReads - nGoodReads) << endl;

	if(nReads > 0){
		cout << "Remaining reads                   " << nGoodReads << "   (" << fixed <<
		        setprecision(2) << 100 * nGoodReads / nReads << "% of input reads)\n" << endl;
	}
}


void printCompletedMessage(const Options &o){
	
	using namespace std;
	using namespace flexbar;
	
	cout << "Flexbar completed ";
	
	if(o.barDetect != BOFF)                      cout << "barcode";
	if(o.barDetect == WITHIN_READ_REMOVAL)       cout << " removal within reads";
	if(o.barDetect == WITHIN_READ)               cout << " detection within reads";
	if(o.barDetect == BARCODE_READ)              cout << " detection with separate reads";
	if(o.barDetect != BOFF && o.adapRm != AOFF)  cout << " and ";
	if(o.barDetect == BOFF && o.adapRm == AOFF)  cout << "basic processing";
	if(o.adapRm != AOFF)                         cout << "adapter removal";
	cout << ".\n" << endl;
}


// #include <seqan/find.h>
// #include <seqan/align.h>

void performTest(){
	
	using seqan::CharString;
	
	// CharString haystack = "ATGGATTGCG";
	// CharString needle   = "ATGCAT";
	// 
	// seqan::Finder<CharString> finder(haystack);
	// seqan::Pattern<CharString, seqan::DPSearch<seqan::SimpleScore, seqan::FindInfix> > pattern(needle, seqan::SimpleScore(0, -1, -7));
	// 
	// while (find(finder, pattern, -2)){
	//     while (findBegin(finder, pattern, getScore(pattern))){
	//         cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << endl;
	//         //cout << end(finder) << endl; //',' << position(pattern) << endl;
	// 	}
	// }
	// 
	// cout << "------" << endl;
	// 
	// clear(finder);
	// seqan::Pattern<CharString, seqan::AbndmAlgo > pattern2(needle, -2);
	// 
	// //seqan::Score<int> sc(0,-3,-2);  // = scoringScheme(pattern2);
	// //setScoringScheme(pattern2, sc);
	// 
	// while (find(finder, pattern2)){
	//     while (findBegin(finder, pattern2, getScore(pattern2))){
	//         cout << '[' << beginPosition(finder) << ',' << endPosition(finder) << ")\t" << infix(finder) << endl;
	// 	}
	// }
	
}


#endif /* FLEXBAR_FLEXBAR_H_ */
