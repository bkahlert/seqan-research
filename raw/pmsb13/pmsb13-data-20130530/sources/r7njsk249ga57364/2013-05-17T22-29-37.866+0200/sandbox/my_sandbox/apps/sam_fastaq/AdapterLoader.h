/*
 *   AdapterLoader.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_ADAPTERLOADER_H_
#define FLEXBAR_ADAPTERLOADER_H_

#include <sstream>
#include <string>

#include <tbb/pipeline.h>
#include <tbb/concurrent_vector.h>

#include <seqan/basic.h>

#include "Enums.h"
#include "SequencingRead.h"
#include "SequenceConverter.h"


typedef std::pair< SequencingRead<seqan::CharString, seqan::CharString>*, tbb::atomic<unsigned long> > TAdapter;


/* This class will store each processed read plus it's ID in a vector. */

template <class TString, class TIDString>
class AdapterLoader : public tbb::filter{

private:
	
	flexbar::FileFormat m_fformat;
	tbb::concurrent_vector<TAdapter> adapters;
	
public:
	
	AdapterLoader(flexbar::FileFormat fformat) : tbb::filter(serial){
		m_fformat = fformat;
	};
	
	virtual ~AdapterLoader(){};
	
	
	void* operator()( void* item ){
		
		SequencingRead<TString, TIDString> *myRead = static_cast< SequencingRead<TString, TIDString>* >(item);
		
		if(m_fformat == flexbar::CSFASTA || m_fformat == flexbar::CSFASTQ){
			TString csRead = SequenceConverter<TString>::getInstance()->basepairSpaceToColorSpace(myRead->getSequence());
			myRead->setSequence(csRead);
		}
		
		TAdapter adap;
		adap.first = myRead;
		adapters.push_back(adap);
		
		return NULL;
	};
	
	
	tbb::concurrent_vector<TAdapter> getAdapters(){
		return adapters;
	}
	
	
	void setAdapters(tbb::concurrent_vector<TAdapter> &adapterVec){
		adapters = adapterVec;
	}
	
	
	void printAdapters(std::string adapterName) const {
		using namespace std;
		
		const unsigned int maxSpaceLen = 23;
		
		cout << adapterName << ":" << string(maxSpaceLen - 8, ' ') << "Sequence:" << "\n";
		
		for(unsigned int i=0; i < adapters.size(); ++i){
			seqan::CharString seqTag = adapters.at(i).first->getSequenceTag();
			
			int whiteSpaceLen = maxSpaceLen - length(seqTag);
			if(whiteSpaceLen < 2) whiteSpaceLen = 2;
			
			string whiteSpace = string(whiteSpaceLen, ' ');
			
			cout << seqTag << whiteSpace << adapters.at(i).first->getSequence() << "\n";
		}
		cout << endl;
	}
	
};

#endif /* FLEXBAR_ADAPTERLOADER_H_ */
