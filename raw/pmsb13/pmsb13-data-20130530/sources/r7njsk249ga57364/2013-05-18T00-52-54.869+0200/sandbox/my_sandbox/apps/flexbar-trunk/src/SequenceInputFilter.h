/*
 *   SequenceInputFilter.h
 *
 *   Authors: mat and jtr
 */

#ifndef FLEXBAR_SEQUENCEINPUTFILTER_H_
#define FLEXBAR_SEQUENCEINPUTFILTER_H_

#include <string>
#include <iostream>
#include <stdexcept>

#include <tbb/pipeline.h>

// #include <seqan/sequence.h>
// #include <seqan/seq_io.h>

#include "Options.h"
#include "Enums.h"
#include "SequencingRead.h"


// This class reads a (CS)FASTA/Q file and builds instances for each SequencingRead.

template <class TString, class TIDString>
class SequenceInputFilter : public tbb::filter {

private:
    
	std::ifstream fstrm;
	
	const flexbar::QualityType m_qualType;
	flexbar::FileFormat m_format;
	TIDString m_nextTag;
	
	const bool m_switch2Fasta, m_preProcess;
	int m_maxUncalled, m_preTrimBegin, m_preTrimEnd, m_prePhredTrim;
	tbb::atomic<unsigned long> m_nrReads, m_nLowPhred;
	
public:
	
	SequenceInputFilter(const Options &o, const std::string filePath, const bool fastaFormat, const bool preProcess) :
		
		filter(serial_in_order),
		m_preProcess(preProcess),
		m_switch2Fasta(o.switch2Fasta),
		m_qualType(o.qual){
		
		m_format       = o.fformat;
		m_maxUncalled  = o.maxUncalled;
		m_preTrimBegin = o.cutLen_begin;
		m_preTrimEnd   = o.cutLen_end;
		
		m_prePhredTrim = 0;
		m_nrReads      = 0;
		m_nLowPhred    = 0;
		
		m_nextTag = "";
		
		using namespace std;
		using namespace flexbar;
		
		if(fastaFormat){
			m_format = FASTA;
		}
		else if(m_switch2Fasta){
			if(m_format == FASTA)   m_format = FASTQ;
			if(m_format == CSFASTA) m_format = CSFASTQ;
		}
		
		// seqan::SequenceStream seqStream("myTest.fa.gz");
		// seqan::CharString id, seq;
		// 
		// readRecord(id, seq, seqStream);
		// cout << id << ": " << seq << endl;
		// int res = readRecord(id, seq, qual, seqStream);
		
		fstrm.open(filePath.c_str(), ios_base::in);
		
		if(! fstrm.is_open()){
			cerr << "Error opening file: " << filePath << endl;
			exit(1);
		}
	};
	
	
	virtual ~SequenceInputFilter(){
		fstrm.close();
	};
	
	
	void setPrePhredTrim(const int nr, const bool showText){
		
		using namespace std;
		using namespace flexbar;
		
		if(nr > 0) m_prePhredTrim = nr;
		else       return;
		
		if(showText) cout << "Trimming reads from 3' end to phred quality ";
		
		switch(m_qualType){
			case SANGER:{
				m_prePhredTrim += 33;
				
				if(showText) cout << nr << " (" << m_prePhredTrim << ")" << endl;
				break;
			}
			case SOLEXA:{
				m_prePhredTrim += 59;
				
				if(showText) cout << nr << " (" << m_prePhredTrim << ")" << endl;
				break;
			}
			case ILLUMINA13:{
				m_prePhredTrim += 64;
				
				if(showText) cout << nr << " (" << m_prePhredTrim << ")" << endl;
				break;
			}
		}
	}
	
	
	unsigned long getNrLowPhredReads() const {
		return m_nLowPhred;
	}
	
	
	unsigned long getNrProcessedReads() const {
		return m_nrReads;
	}
	
	
	std::string readLine(){
		std::string text;
		
		if(fstrm.good())
		getline(fstrm, text);
		
		return text;
	}
	
	
	// Core method for reading and parsing FASTA/FASTQ input files.
	// @return: single SequencingRead<TString, TIDString> or NULL if no more reads in file or error.
	
	void* getRead(bool &isUncalled){
		
		using namespace std;
		using namespace flexbar;
		
		using seqan::prefix;
		using seqan::suffix;
		using seqan::length;
		
		SequencingRead<TString, TIDString> *myRead = NULL;
		
		TString source = "", quality = "", dummy = "";
		TIDString tag = "";
		
		if(! fstrm.eof()){
			
			isUncalled = false;
			
			try{
				// FastA parsing
				if(m_format == FASTA || m_format == CSFASTA){
					
					// tag line is read in previous iteration
					if(m_nextTag == "") tag = readLine();
					else                tag = m_nextTag;
					
					if(length(tag) > 0){
						if(seqan::isNotEqual(getValue(tag, 0), '>')){
							stringstream error;
							error << "Incorrect FASTA entry, missing > on new line. Input: " << tag << endl;
							throw runtime_error(error.str());
						}
						else tag = suffix(tag, 1);
						
						if(length(tag) == 0){
							stringstream error;
							error << "Incorrect FASTA entry, missing read name after > symbol." << endl;
							throw runtime_error(error.str());
						}
					}
					else return NULL;
					
					
					source = readLine();
					
					if(length(source) < 1){
						stringstream error;
						error << "Empty FASTA entry, found tag without read! Tag: " << tag << endl;
						throw runtime_error(error.str());
					}
					
					
					m_nextTag = readLine();
					
					// fasta files with sequences splitted over several lines
					while(fstrm.good() && length(m_nextTag) > 0 && seqan::isNotEqual(getValue(m_nextTag, 0), '>')){
						append(source, m_nextTag);
						m_nextTag = readLine();
					}
					
					if(m_preProcess){
						isUncalled = isUncalledSequence(source);
						
						if(m_preTrimBegin > 1 && length(source) > m_preTrimBegin){
							source = suffix(source, m_preTrimBegin);

							if(m_format == CSFASTA) insert(source, 0, "T");
						}
						
						if(m_preTrimEnd > 1 && length(source) > m_preTrimEnd){
							source = prefix(source, length(source) - m_preTrimEnd);
						}
					}
					
					myRead = new SequencingRead<TString, TIDString>(source, tag);
					
					++m_nrReads;
				}
				
				// FastQ parsing
				else{
					
					source = readLine();
					
					if(length(source) > 0){
						if(seqan::isNotEqual(getValue(source, 0), '@')){
							stringstream error;
							error << "Incorrect FASTQ entry, missing @ on new line. Input: " << source << endl;
							throw runtime_error(error.str());
						}
						else tag = suffix(source, 1);
						
						if(length(tag) == 0){
							stringstream error;
							error << "Incorrect FASTQ entry, missing read name after @ symbol." << endl;
							throw runtime_error(error.str());
						}
					}
					else return NULL;
					
					source = readLine();
					
					if(length(source) < 1){
						stringstream error;
						error << "Empty FASTQ entry, found tag without read! Tag: " << tag << endl;
						throw runtime_error(error.str());
					}
					
					
					dummy = readLine();
					if(length(dummy) == 0 || seqan::isNotEqual(getValue(dummy, 0), '+')){
							stringstream error;
							error << "Incorrect FASTQ entry, missing + line. Tag: " << tag << endl;
							throw runtime_error(error.str());
					}
					
					quality = readLine();
					
					// in case CSFASTQ format has same quality and read length it will be trimmed
					if(m_format == CSFASTQ){
						if(length(quality) == length(source)){
							quality = suffix(quality, 1);
						}
					}
					
					if(length(quality) < 1){
						stringstream error;
						error << "Empty FASTQ entry, found read without quality values! Tag: " << tag << endl;
						throw runtime_error(error.str());
					}
					
					
					if(m_preProcess){
						isUncalled = isUncalledSequence(source);
						
						if(m_preTrimBegin > 1 && length(source) > m_preTrimBegin){
							source = suffix(source, m_preTrimBegin);

							if(m_format == CSFASTQ){
								insert(source, 0, "T");

								quality = suffix(quality, m_preTrimBegin - 1);
							}
							else quality = suffix(quality, m_preTrimBegin);
						}
						
						if(m_preTrimEnd > 1 && length(source) > m_preTrimEnd){
							source = prefix(source, length(source) - m_preTrimEnd);

							if(m_format == CSFASTQ){
								quality = prefix(quality, (length(quality) - m_preTrimEnd) - 1);
							}
							else quality = prefix(quality, length(quality) - m_preTrimEnd);
						}
						
						// filtering based on phred quality
						if(m_prePhredTrim != 0){
							typename seqan::Iterator<TString >::Type it    = seqan::begin(quality);
							typename seqan::Iterator<TString >::Type itEnd = seqan::end(quality);
						
							--itEnd;
						
							unsigned int n = length(quality);
						
							bool nChanged = false;
						
							while(itEnd != it){
								if(static_cast<int>(*itEnd) >= m_prePhredTrim) break;
								--n;
								--itEnd;
							
								if(! nChanged){
									m_nLowPhred++;
									nChanged = true;
								}
							}
							source = prefix(source, n);
						
							if(m_format == CSFASTQ) --n;
							quality = prefix(quality, n);
						}
					}
					
					if(m_switch2Fasta) myRead = new SequencingRead<TString, TIDString>(source, tag);
					else               myRead = new SequencingRead<TString, TIDString>(source, tag, quality);
					
					++m_nrReads;
				}
				
				return myRead;
			}
			catch(exception &e){
				cerr << "\n\n" << e.what() << "\nProgram execution aborted.\n" << endl;
				fstrm.close();
				
				exit(1);
			}
		}
		
		// end of stream
		else return NULL;
	}
	
	
	// returns TRUE if read contains too many uncalled bases
	bool isUncalledSequence(TString source){
		int n = 0;
		
		typename seqan::Iterator<TString >::Type it, itEnd;
		
		it    = seqan::begin(source);
		itEnd = seqan::end(source);
		
		while(it != itEnd){
			 if(*it == '.' || *it == 'N') n++;
			 ++it;
		}
		
		return(n > m_maxUncalled);
	}
 	
	
	// override
	void* operator()(void*){
		
		bool isUncalled = false;
		return getRead(isUncalled);
	}
	
};

#endif /* FLEXBAR_SEQUENCEINPUTFILTER_H_ */
