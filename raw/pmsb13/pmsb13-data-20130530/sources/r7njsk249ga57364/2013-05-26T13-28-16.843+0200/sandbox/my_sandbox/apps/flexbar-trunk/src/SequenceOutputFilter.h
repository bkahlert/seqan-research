/*
 *   SequenceOutputFilter.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_SEQUENCEOUTPUTFILTER_H_
#define FLEXBAR_SEQUENCEOUTPUTFILTER_H_

#include <fstream>

#include <tbb/concurrent_vector.h>

#include "Enums.h"
#include "SequencingRead.h"


// This class writes sequencing reads in specified format to a file.

template <typename TString, typename TIDString>
class SequenceOutputFilter {

private:
	
	unsigned int m_minLength, m_cutLen_read;
	
	tbb::concurrent_vector<unsigned long> *m_lengthDist;
	
	std::ofstream m_targetStream;
	std::string m_filePath;
	
	tbb::atomic<unsigned long> m_countGood;
	flexbar::FileFormat m_format;
	
	bool m_writeLenDist;
	
public:
	
	SequenceOutputFilter(const std::string& filePath, Options &o){
		using namespace std;
		
		m_format       = o.format;
		m_minLength    = o.min_readLen;
		m_cutLen_read  = o.cutLen_read;
		m_writeLenDist = o.writeLengthDist;
		
		m_filePath = filePath;
		
		m_countGood = 0;
		
		m_lengthDist = new tbb::concurrent_vector<unsigned long>(flexbar::MAX_READLENGTH + 1, 0);
		
		m_targetStream.open(m_filePath.c_str(), ios_base::out | ios_base::binary);
		
		if(! m_targetStream.is_open() || ! m_targetStream.good()){
			cerr << "Error opening file: " << m_filePath << "\n";
			exit(1);
		}
	};
	
	
	virtual ~SequenceOutputFilter(){
		m_targetStream.close();
		delete m_lengthDist;
	};
	
	
	const std::string& getFileName() const{
		return m_filePath;
	}
	
	
	void writeLengthDist() const {
		using namespace std;
		
		string fname = m_filePath + ".lengthdist";
		fstream lstream;
		
		lstream.open(fname.c_str(), ios_base::out | ios_base::binary);
		
		if(! lstream.is_open()){
			cerr << "Error opening File: " << fname << "\n";
		}
		else{
			lstream << "Readlength\tCount" << "\n";
			
			for (int i = 0; i <= flexbar::MAX_READLENGTH; ++i){
				if(m_lengthDist->at(i) > 0)
					lstream << i << "\t" << m_lengthDist->at(i) << "\n";
			}
			lstream.close();
		}
	}
	
	
	void writeFastString(const SequencingRead<TString, TIDString>& myRead){
		
		using namespace std;
		using namespace flexbar;
		
		switch(m_format){
			case FASTQ: 
			case CSFASTQ:
			m_targetStream << "@" << myRead.getSequenceTag() << "\n" << myRead.getSequence() << "\n+\n" << myRead.getQuality() << "\n";
			break;
			
			case FASTA:
			case CSFASTA:
			m_targetStream << ">" << myRead.getSequenceTag() << "\n" << myRead.getSequence() << "\n";
			break;
		}
	}
	
	
	unsigned long getNrGoodReads() const {
		return m_countGood;
	}
	
	
	void *writeRead(void *item){
		
		using namespace std;
		using namespace flexbar;
		
		if(item){
			SequencingRead<TString, TIDString> *myRead = static_cast< SequencingRead<TString, TIDString>* >(item);
			
			unsigned int readLength = length(myRead->getSequence());
			
			if(m_cutLen_read > 1 && m_cutLen_read >= m_minLength && m_cutLen_read < readLength){
				
				myRead->setSequence(prefix(myRead->getSequence(), m_cutLen_read));
				
				if(m_format == FASTQ){
					myRead->setQuality(prefix(myRead->getQuality(), m_cutLen_read));
				}
				else if(m_format == CSFASTQ){
					myRead->setQuality(prefix(myRead->getQuality(), m_cutLen_read - 1));
				}
			}
			
			++m_countGood;
			
			// store read length distribution
			if(m_writeLenDist && readLength <= MAX_READLENGTH)
				m_lengthDist->at(readLength)++;
			else if(m_writeLenDist)
				cerr << "\nCompile Flexbar with larger max read length to get correct length dist.\n" << endl;
			
			writeFastString(*myRead);
		}
		
		return NULL;
	}
	
};

#endif /* FLEXBAR_SEQUENCEOUTPUTFILTER_H_ */
