/*
 *   MultiplexedInputFilter.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_MULTIPLEXEDOUTPUTFILTER_H_
#define FLEXBAR_MULTIPLEXEDOUTPUTFILTER_H_

#include <tbb/pipeline.h>
#include <tbb/concurrent_vector.h>

#include <seqan/basic.h>

#include "Options.h"
#include "Enums.h"
#include "MultiplexedRead.h"
#include "SequenceOutputFilter.h"
#include "OutputFileStruct.h"
#include "AdapterLoader.h"


/* This class will process a MultiplexedRead and write it to a file
   depending on the runtype: single-end, paired-end and/or barcoded. */

template <typename TString, typename TIDString>
class MultiplexedOutputFilter : public tbb::filter {

private:
	
	int m_mapsize;
	const int m_minLength, m_cutLen_read;
	
	const bool m_writeSingleReads;
	const std::string m_target;
	const flexbar::FileFormat m_format;
	const flexbar::RunType m_runType;
	
	typedef SequenceOutputFilter<TString, TIDString> TOutputFilter;
	typedef OutputFileStruct<TString, TIDString> filters;
	
	filters *m_outputMap;
	
	tbb::concurrent_vector<TAdapter> *m_adapters, *m_barcodes;
	
public:
	
	MultiplexedOutputFilter(Options &o, tbb::concurrent_vector<TAdapter> *barcodes, tbb::concurrent_vector<TAdapter> *adapters) :
		
		filter(serial_in_order),
		m_target(o.targetName),
		m_format(o.format),
		m_runType(o.runType),
		m_minLength(o.min_readLen),
		m_cutLen_read(o.cutLen_read),
		m_writeSingleReads(o.writeSingleReads){
		
		using namespace std;
		using namespace flexbar;
		
		m_adapters  = adapters;
		m_barcodes  = barcodes;
		
		m_mapsize = 0;
		
		
		switch(m_runType){
			
			case PAIRED_BARCODED:{
				
				m_mapsize = m_barcodes->size() + 1;
				m_outputMap = new filters[m_mapsize];
				
				stringstream ss;
				
				for(unsigned int i = 0; i < m_barcodes->size(); ++i){
					
					ss << m_target << "_barcode_" << m_barcodes->at(i).first->getSequenceTag() << "_1" << toFormatString(m_format);
					TOutputFilter *of1 = new TOutputFilter(ss.str(), o);
					ss.str("");
					ss.clear();
					
					ss << m_target << "_barcode_" << m_barcodes->at(i).first->getSequenceTag() << "_2"<< toFormatString(m_format);
					TOutputFilter *of2 = new TOutputFilter(ss.str(), o);
					ss.str("");
					ss.clear();
					
					filters& f = m_outputMap[i + 1];
					f.f1       = of1;
					f.f2       = of2;
					
					if(m_writeSingleReads){
						ss << m_target << "_barcode_" << m_barcodes->at(i).first->getSequenceTag() << "_1_single" << toFormatString(m_format);
						TOutputFilter *osingle1 = new TOutputFilter(ss.str(), o);
						ss.str("");
						ss.clear();
						
						ss << m_target << "_barcode_" << m_barcodes->at(i).first->getSequenceTag() << "_2_single"<< toFormatString(m_format);
						TOutputFilter *osingle2 = new TOutputFilter(ss.str(), o);
						ss.str("");
						ss.clear();
						
						f.single1 = osingle1;
						f.single2 = osingle2;
					}
				}
				
				ss << m_target << "_barcode_unassigned_1"<< toFormatString(m_format);
				TOutputFilter *of1 = new TOutputFilter(ss.str(), o);
				ss.str("");
				ss.clear();
				
				ss << m_target << "_barcode_unassigned_2"<< toFormatString(m_format);
				TOutputFilter *of2 = new TOutputFilter(ss.str(), o);
				ss.str("");
				ss.clear();
				
				filters& f = m_outputMap[0];
				f.f1       = of1;
				f.f2       = of2;
				
				if(m_writeSingleReads){
					ss << m_target << "_barcode_unassigned_1_single"<< toFormatString(m_format);
					TOutputFilter *osingle1 = new TOutputFilter(ss.str(), o);
					ss.str("");
					ss.clear();
					
					ss << m_target << "_barcode_unassigned_2_single"<< toFormatString(m_format);
					TOutputFilter *osingle2 = new TOutputFilter(ss.str(), o);
					ss.str("");
					ss.clear();
					
					f.single1 = osingle1;
					f.single2 = osingle2;
				}
				break;
			}
			
			case PAIRED:{
				
				m_mapsize = 1;
				m_outputMap = new filters[m_mapsize];
				stringstream ss;
				
				ss << m_target << "_1"<< toFormatString(m_format);
				TOutputFilter *of1 = new TOutputFilter(ss.str(), o);
				ss.str("");
				ss.clear();
				
				ss << m_target << "_2"<< toFormatString(m_format);
				TOutputFilter *of2 = new TOutputFilter(ss.str(), o);
				ss.str("");
				ss.clear();
				
				filters& f = m_outputMap[0];
				f.f1       = of1;
				f.f2       = of2;
				
				if(m_writeSingleReads){
					ss << m_target << "_1_single" << toFormatString(m_format);
					TOutputFilter *osingle1 = new TOutputFilter(ss.str(), o);
					ss.str("");
					ss.clear();
					
					ss << m_target << "_2_single"<< toFormatString(m_format);
					TOutputFilter *osingle2 = new TOutputFilter(ss.str(), o);
					ss.str("");
					ss.clear();
					
					f.single1 = osingle1;
					f.single2 = osingle2;
				}
				break;
			}
			
			case SINGLE:{
				
				m_mapsize = 1;
				
				m_outputMap = new filters[m_mapsize];
				stringstream ss;
				
				ss << m_target << toFormatString(m_format);
				TOutputFilter *of1 = new TOutputFilter(ss.str(), o);
				ss.str("");
				ss.clear();
				
				filters& f = m_outputMap[0];
				f.f1 = of1;
				
				break;
			}
			
			case SINGLE_BARCODED:{
				
				m_mapsize = m_barcodes->size() + 1;
				m_outputMap = new filters[m_mapsize];
				
				stringstream ss;
				
				for(unsigned int i=0; i < m_barcodes->size(); ++i){
					
					ss << m_target << "_barcode_" << m_barcodes->at(i).first->getSequenceTag() << toFormatString(m_format);
					TOutputFilter *of1 = new TOutputFilter(ss.str(), o);
					ss.str("");
					ss.clear();
					
					filters& f = m_outputMap[i+1];
					f.f1 = of1;
				}
				
				ss << m_target << "_barcode_unassigned" << toFormatString(m_format);
				TOutputFilter *of1 = new TOutputFilter(ss.str(), o);
				ss.str("");
				ss.clear();
				
				filters& f = m_outputMap[0];
				f.f1 = of1;
				
				break;
			}
		}
	}
	
	
	virtual ~MultiplexedOutputFilter(){
		delete[] m_outputMap;
	};
	
	
	void* operator()(void* item) {
		
		using namespace flexbar;
		
		MultiplexedRead<TString, TIDString> *myRead = static_cast< MultiplexedRead<TString, TIDString>* >(item);
		
		bool l1ok = false, l2ok = false;
		
		switch(m_runType){
			
			case SINGLE:
			case SINGLE_BARCODED:{
				
				if(myRead->m_r1 != NULL){
					if(length(myRead->m_r1->getSequence()) >= m_minLength){
						
						m_outputMap[myRead->m_barcode_id].f1->writeRead(myRead->m_r1);
					}
					else m_outputMap[myRead->m_barcode_id].m_nShort_1++;
				}
				break;
			}
			
			case PAIRED:
			case PAIRED_BARCODED:{
				
				if(myRead->m_r1 != NULL && myRead->m_r2 != NULL){
					
					// now check if both reads have minLength
					if(length(myRead->m_r1->getSequence()) >= m_minLength) l1ok = true;
					if(length(myRead->m_r2->getSequence()) >= m_minLength) l2ok = true;
					
					if(l1ok && l2ok){
						m_outputMap[myRead->m_barcode_id].f1->writeRead(myRead->m_r1);
						m_outputMap[myRead->m_barcode_id].f2->writeRead(myRead->m_r2);
					}
					else if(m_writeSingleReads && l1ok && ! l2ok){
						m_outputMap[myRead->m_barcode_id].single1->writeRead(myRead->m_r1);
					}
					else if(m_writeSingleReads && ! l1ok && l2ok){
						m_outputMap[myRead->m_barcode_id].single2->writeRead(myRead->m_r2);
					}
					
					if(! l1ok) m_outputMap[myRead->m_barcode_id].m_nShort_1++;
					if(! l2ok) m_outputMap[myRead->m_barcode_id].m_nShort_2++;
				}
			}
		}
		
		delete myRead;
		
		return NULL;
	}
	
	
	void writeLengthDist(){
		for(unsigned int i = 0; i < m_mapsize; i++){
			m_outputMap[i].f1->writeLengthDist();
			if(m_outputMap[i].f2 != NULL) m_outputMap[i].f2->writeLengthDist();
		}
	}
	
	
	unsigned long getNrGoodReads(){
		unsigned long good = 0;
		
		for(unsigned int i = 0; i < m_mapsize; i++){
			good += m_outputMap[i].f1->getNrGoodReads();
			if(m_outputMap[i].f2 != NULL) good += m_outputMap[i].f2->getNrGoodReads();
		}
		return good;
	}
	
	
	unsigned long getNrShortReads(){
		unsigned long nShortReads = 0;
		
		for(unsigned int i = 0; i < m_mapsize; i++){
			nShortReads += m_outputMap[i].m_nShort_1;
			
			if(m_runType == flexbar::PAIRED_BARCODED || m_runType == flexbar::PAIRED){
				nShortReads += m_outputMap[i].m_nShort_2;
			}
		}
		return nShortReads;
	}
	
	
	void printAdapterRemovalStats(){
		using namespace std;
		
		cout << "Adapter removal statistics\n";
		cout << "==========================\n";
		
		const unsigned int maxSpaceLen = 18;
		
		cout << "Adapter:" << string(maxSpaceLen - 8, ' ') << "Removed:" << "\n";
		
		for(unsigned int i = 0; i < m_adapters->size(); i++){
			seqan::CharString seqTag = m_adapters->at(i).first->getSequenceTag();
			
			int whiteSpaceLen = maxSpaceLen - length(seqTag);
			if(whiteSpaceLen < 2) whiteSpaceLen = 2;
			
			string whiteSpace = string(whiteSpaceLen, ' ');
			
			cout << seqTag << whiteSpace << m_adapters->at(i).second << "\n";
		}
		cout << endl;
	}
	
	
	void printFileSummary(){
		
		using namespace std;
		
		cout << "Output file statistics\n";
		cout << "======================\n";
		
		for(unsigned int i = 0; i < m_mapsize; i++){
			cout << "Read file:               " << m_outputMap[i].f1->getFileName()    << "\n";
			cout << "  written reads          " << m_outputMap[i].f1->getNrGoodReads() << "\n";
			cout << "  skipped short reads    " << m_outputMap[i].m_nShort_1           << endl;
			
			if(m_runType == flexbar::PAIRED_BARCODED || m_runType == flexbar::PAIRED){
				cout << "Read file2:              " << m_outputMap[i].f2->getFileName()    << "\n";
				cout << "  written reads          " << m_outputMap[i].f2->getNrGoodReads() << "\n";
				cout << "  too short reads        " << m_outputMap[i].m_nShort_2           << endl;
			}
		}
		cout << endl;
	}
	
};

#endif /* FLEXBAR_MULTIPLEXEDOUTPUTFILTER_H_ */
