/*
 *   SequencingRead.h
 *
 *   Author: mat and jtr
 */

#ifndef FLEXBAR_SEQUENCINGREAD_H_
#define FLEXBAR_SEQUENCINGREAD_H_

#include <seqan/basic.h>


/* A Sequencing read consists of a nucleotide sequence (color or basepair space),
   a sequence name and optionally a quality string plus the quality scaling. */

template <class TString, class TIDString>
class SequencingRead {

private:
	
	TString m_seq;
	TIDString m_tag, m_qual;
	
public:
	
	SequencingRead()
    	: m_tag(),
      	  m_seq(){
	}
	
	SequencingRead(const TString& source, const TIDString& sequence_tag)
		: m_tag(sequence_tag),
		  m_seq(source){
	}
	
	SequencingRead(const TString& source, const TIDString& sequence_tag, const TIDString& qual)
		: m_tag(sequence_tag),
		  m_seq(source),
		  m_qual(qual){
	}
	
	
	void setSequenceTag(const TString& tag){
		m_tag = tag;
	}
	
	void setSequence(const TString& seq){
		m_seq = seq;
	}
	
	void setQuality(const TString& qual){
		m_qual = qual;
	}
	
	
	const TIDString& getSequenceTag() const {
		return m_tag;
	}
	
	const TString& getSequence() const {
		return m_seq;
	}
	
	const TIDString& getQuality() const{
		return m_qual;
	}
	
	
	virtual ~SequencingRead(){};
	
};

#endif /* FLEXBAR_SEQUENCINGREAD_H_ */
