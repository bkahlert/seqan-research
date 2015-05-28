/*
 *   MultiplexedInputFilter.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_MULTIPLEXEDINPUTFILTER_H_
#define FLEXBAR_MULTIPLEXEDINPUTFILTER_H_

#include <tbb/pipeline.h>

#include "Options.h"
#include "MultiplexedRead.h"
#include "SequenceInputFilter.h"


/* This class handles up to 3 file sources (paired and barcode reads)
   and creates a MultiplexedRead depending on the run type. */

template <typename TString, typename TIDString>
class MultiplexedInputFilter : public tbb::filter {

private:
	
	const bool m_isPaired, m_useBarcodeRead, m_useNumberTag;
	tbb::atomic<unsigned long> m_uncalled, m_uncalledPairs, m_tagCounter;
	
	SequenceInputFilter<TString, TString > *m_f1, *m_f2, *m_b;
	
public:
	
	MultiplexedInputFilter(const Options &o) :
		
		filter(serial_in_order),
		m_useNumberTag(o.useNumberTag),
		m_isPaired(o.runType == flexbar::PAIRED || o.runType == flexbar::PAIRED_BARCODED),
		m_useBarcodeRead(o.barDetect == flexbar::BARCODE_READ){
		
		m_uncalled      = 0;
		m_uncalledPairs = 0;
		m_tagCounter    = 0;
		
		m_f1 = new SequenceInputFilter<TString, TString>(o, o.readsFile, false, true);
		m_f1->setPrePhredTrim(o.phred_preQual, true);
		
		m_f2 = NULL;
		m_b  = NULL;
		
		if(m_isPaired){
			m_f2 = new SequenceInputFilter<TString, TString>(o, o.readsFile2, false, true);
			m_f2->setPrePhredTrim(o.phred_preQual, false);
		}
		
		if(m_useBarcodeRead){
			m_b = new SequenceInputFilter<TString, TString>(o, o.barReadsFile, false, false);
		}
	}
	
	
	virtual ~MultiplexedInputFilter(){
		delete m_f1;
		delete m_f2;
		delete m_b;
	}
	
	
	void* operator()(void*){
		
		using namespace std;
		
		SequencingRead<TString, TIDString> *myRead1 = NULL, *myRead2 = NULL, *myBarcode = NULL;
		
		bool uncalled = true, uncalled2 = true, uBR = true;
		
		// single read input
		if(! m_isPaired){
			
			while(uncalled){
				myRead1 = static_cast< SequencingRead<TString, TIDString>* >(m_f1->getRead(uncalled));
				
				if(m_useBarcodeRead) myBarcode = static_cast< SequencingRead<TString, TIDString>* >(m_b->getRead(uBR));
				
				if(myRead1 == NULL) return NULL;
				
				else if(m_useBarcodeRead && myBarcode == NULL){
					cerr << "Error: read without barcode read, or file reading error!\n" << endl;
					exit(1);
				}
				
				if(uncalled){
					++m_uncalled;
					delete myRead1;
					delete myBarcode;
				}
			}
		}
		
		// paired read input
		else{
			
			// find paired read without uncalled bases
			while(uncalled || uncalled2){
				
				myRead1 = static_cast< SequencingRead<TString, TIDString>* >(m_f1->getRead(uncalled));
				myRead2 = static_cast< SequencingRead<TString, TIDString>* >(m_f2->getRead(uncalled2));
				
				if(m_useBarcodeRead) myBarcode = static_cast< SequencingRead<TString, TIDString>* >(m_b->getRead(uBR));
				
				// end of files reached
				if(myRead1 == NULL && myRead2 == NULL) return NULL;
				
				else if(myRead1 == NULL || myRead2 == NULL){
					cerr << "Error: single read in paired mode, or file reading error!\n" << endl;
					exit(1);
				}
				else if(m_useBarcodeRead && myBarcode == NULL){
					cerr << "Error: reads without barcode read or file reading error!\n" << endl;
					exit(1);
				}
				
				if(uncalled || uncalled2){
					++m_uncalledPairs;
					if(uncalled)  ++m_uncalled;
					if(uncalled2) ++m_uncalled;
					
					delete myRead1;
					delete myRead2;
					delete myBarcode;
				}
			}
		}
		
		if(m_useNumberTag){
			stringstream converter;
			converter << ++m_tagCounter;
			TString tagCount = converter.str();
			
			myRead1->setSequenceTag(tagCount);
			if(m_isPaired) myRead2->setSequenceTag(tagCount);
			if(m_useBarcodeRead) myBarcode->setSequenceTag(tagCount);
		}
		
		return new MultiplexedRead<TString, TIDString>(myRead1, myRead2, myBarcode);
	}
	
	
	unsigned long getNrUncalledReads() const{
		return m_uncalled;
	}
	
	
	unsigned long getNrUncalledPairedReads() const{
		return m_uncalledPairs;
	}
	
	
	unsigned long getNrProcessedReads() const{
		if(m_isPaired) return m_f1->getNrProcessedReads() + m_f2->getNrProcessedReads();
		else           return m_f1->getNrProcessedReads();
	}
	
	
	unsigned long getNrLowPhredReads() const {
		if(m_isPaired) return m_f1->getNrLowPhredReads() + m_f2->getNrLowPhredReads();
		else           return m_f1->getNrLowPhredReads();
	}
	
};

#endif /* FLEXBAR_MULTIPLEXEDINPUTFILTER_H_ */
