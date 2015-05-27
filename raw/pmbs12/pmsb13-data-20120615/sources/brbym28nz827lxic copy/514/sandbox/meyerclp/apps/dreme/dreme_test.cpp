#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>

#include "dreme.h"

SEQAN_DEFINE_TEST(test_ChIPSeq_fasta_print){

	using namespace seqan;
	Seq sequences;

	char buffer[1023];
    strcpy(buffer, SEQAN_PATH_TO_ROOT());
    strcat(buffer, "/sandbox/meyerclp/tests/ChIPSeq/example.fasta");
    
	readFastA(sequences, buffer);

	SEQAN_ASSERT_EQ(length(sequences.ids), length(sequences.seqs));

    SEQAN_ASSERT_EQ(sequences.ids[0], "first sequence");
    SEQAN_ASSERT_EQ(sequences.ids[1], "second sequence");
    SEQAN_ASSERT_EQ(sequences.ids[2], "third sequence");
    SEQAN_ASSERT_EQ(sequences.ids[3], "fourth sequence");
    SEQAN_ASSERT_EQ(sequences.ids[4], "fifth sequence");
    SEQAN_ASSERT_EQ(sequences.ids[5], "sixth sequence");

    SEQAN_ASSERT_EQ(sequences.seqs[0], "AAATTTTTTTTT");
    SEQAN_ASSERT_EQ(sequences.seqs[1], "CCCGGGA");
    SEQAN_ASSERT_EQ(sequences.seqs[2], "");
    SEQAN_ASSERT_EQ(sequences.seqs[3], "AAAAAA");
    SEQAN_ASSERT_EQ(sequences.seqs[4], "CGA");
    SEQAN_ASSERT_EQ(sequences.seqs[5], "TTTT");
}

SEQAN_DEFINE_TEST(test_ChIPSeq_suffixArray){
	using namespace seqan;

	Seq seq;
	seq.seed=100;
	StringSet<String<Dna5> > test;
	appendValue(test,"ACGTACGT");
	appendValue(test,"TTTTTxxx");
	seq.seqs=test;


    StringSet<String<Dna5> > expected;
    appendValue(expected, "ACGTACGT");
    appendValue(expected, "ACGT");
    appendValue(expected, "CGTACGT");
    appendValue(expected, "CGT");
    appendValue(expected, "GTACGT");
    appendValue(expected, "GT");
    appendValue(expected, "TACGT");
    appendValue(expected, "TTTTTNNN");
    appendValue(expected, "TTTTNNN");
    appendValue(expected, "TTTT");
    appendValue(expected, "TTTNNN");
    appendValue(expected, "TTT");
    appendValue(expected, "TTNNN");
    appendValue(expected, "TT");
    appendValue(expected, "TNNN");
    appendValue(expected, "T");
    appendValue(expected, "NNN");
    appendValue(expected, "NN");
    appendValue(expected, "N");

	initST(seq);


	typedef Index<StringSet<String<Dna5> > > TMyIndex;//Dna5Q
	Iterator<TMyIndex, BottomUp<> >::Type myIterator(seq.SArray);
	for(;!atEnd(myIterator);++myIterator){
		std::cout<<representative(myIterator)<<std::endl;
		
	}
	clear(seq.SArray);
	clear(seq.seqs);

}


SEQAN_DEFINE_TEST(test_ChIPSeq_initExactKmer){
	using namespace seqan;

	Seq seq;
	Seq back;
	seq.seed=100;
	
	StringSet<String<Dna5> > test;
	appendValue(test,"ACGTACGT");
	appendValue(test,"TTTTTTTT");
	appendValue(test,"ACCTACGTTT");
	appendValue(test,"");

	StringSet<String<Dna5> > test2;
	appendValue(test2,"ACG");
	appendValue(test2,"ACGT");
	appendValue(test2,"TTT");

	seq.seqs=test;
	back.seqs=test2;
	seq.c=1;
	back.c=1;
	seq.SeqsNumber=length(seq.seqs);
	back.SeqsNumber=length(back.seqs);
	SEQAN_ASSERT_EQ(seq.SeqsNumber,4);

	SEQAN_ASSERT_EQ(back.SeqsNumber, 3);


	/*typedef std::map<String<Dna5>,std::vector<int> > DnaCounterMap;
	DnaCounterMap seqCount;
	DnaCounterMap backCount;*/
	std::map<String<Dna5>,unsigned int > testMap;
	testMap["ACG"]=2;
	testMap["CGT"]=2;
	testMap["GTA"]=1;
	testMap["TAC"]=2;
	testMap["TTT"]=2;
	testMap["ACC"]=1;
	testMap["CCT"]=1;
	testMap["CTA"]=1;
	testMap["GTT"]=1;

	initST(seq);
	initST(back);
	initExactKmer(seq,back,3,3);

	std::map<String<Dna5>,unsigned int >::iterator MapIterator;
	std::map<String<Dna5>,unsigned int >::iterator MapIterator2;
	MapIterator=seq.seqCounter.begin();
	MapIterator2=testMap.begin();
	SEQAN_ASSERT_EQ(length(seq.seqCounter), length(testMap));
	
	for(;MapIterator != seq.seqCounter.end();++MapIterator,++MapIterator2){
		
		SEQAN_ASSERT_EQ((*MapIterator).first,(*MapIterator2).first);
		SEQAN_ASSERT_EQ((*MapIterator).second,(*MapIterator2).second);
	}


	clear(seq.seqs);
	clear(back.seqs);
	clear(seq.SArray);
	clear(back.SArray);
	
}

SEQAN_DEFINE_TEST(test_ChIPSeq_fact){
	using namespace seqan;
	logFactorial(22);
	SEQAN_ASSERT_EQ(int(exp(logFac[5])+0.5),120);
	SEQAN_ASSERT_EQ(int(exp(logFac[0])+0.5),1);
	SEQAN_ASSERT_EQ(int(exp(logFac[1])+0.5),1);
	SEQAN_ASSERT_EQ(int(exp(logFac[3])+0.5),6);
	SEQAN_ASSERT_EQ(int(exp(logFac[10])+0.5),3628800);

	SEQAN_ASSERT_EQ(calcFET(0,0,0,0),1);

	SEQAN_ASSERT_EQ(double(int((calcFET(3,1,1,1)+0.005)*100))/100,0.53);//auf zwei nachkommastellen runden
	SEQAN_ASSERT_EQ(double(int((calcFET(1,1,1,1)+0.005)*100))/100,0.67);
	SEQAN_ASSERT_EQ(double(int((calcFET(10,0,11,1)+0.005)*100))/100,0.55);

}

SEQAN_DEFINE_TEST(test_ChIPSeq_modifyFET){
	using namespace seqan;
	unsigned int a,b,c,d;
	double pValue=0;

	a=0,b=0,c=0,d=0;
	modifyFET(a,b,c,d,pValue);
	SEQAN_ASSERT_EQ(double(int((pValue+0.005)*100))/100,1.00);

	a=1,b=1,c=1,d=1;
	modifyFET(a,b,c,d,pValue);
	SEQAN_ASSERT_EQ(double(int((pValue+0.005)*100))/100,0.83);

	a=3;
	b=2,c=2,d=2;
	modifyFET(a,b,c,d,pValue);
	
	SEQAN_ASSERT_EQ(double(int((pValue+0.005)*100))/100,0.64);

	a=4;
	b=3;
	c=2;
	d=1;
	modifyFET(a,b,c,d,pValue);
	SEQAN_ASSERT_EQ(double(int((pValue+0.005)*100))/100,0.83);
}


SEQAN_DEFINE_TEST(test_ChIPSeq_generalize){

	using namespace seqan;
	Seq seq;
	Seq back;
	seq.seed=100;
	seq.c=1;
	back.c=1;
	
	
	StringSet<String<Dna5> > test;
	appendValue(test,"ACGTACGT");
	appendValue(test,"TTTTTTTT");
	appendValue(test,"ACCTACGTTT");
	appendValue(test,"");

	StringSet<String<Dna5> > test2;
	appendValue(test2,"ACG");
	appendValue(test2,"ACGT");
	appendValue(test2,"TTT");

	seq.seqs=test;
	back.seqs=test2;

	/*typedef std::map<String<Dna5>,std::vector<int> > DnaCounterMap;
	DnaCounterMap seqCount;
	DnaCounterMap backCount;*/

	seq.SeqsNumber=length(seq.seqs);
	back.SeqsNumber=length(back.seqs);

	initST(seq);
	initST(back);
	initExactKmer(seq,back,3,3);

	//std::map<unsigned int,char> IupacMap;
	//std::map<char,unsigned int> IupacMapReversed;
	//std::map<char,String<Iupac> > IupacMapReplace; 
	//std::map<char,String<Dna5> > IupacMapReplaceReversed;
	//MapIupac(IupacMap, IupacMapReversed, IupacMapReplace,IupacMapReplaceReversed);//IupacMap for generalization

	//std::map<String<Iupac>,unsigned int> generalizedKmer;//unsigned int = estimated counter
	//std::map<String<Iupac>,unsigned int> generalizedKmerBackground;
	//std::map<String<Iupac>,double > generalizedSortedPValueReversed;
	IupacMaps IMaps;
	MapIupac(IMaps);
	GeneralizeKmer("ACG",IMaps,seq,back);
	PrintMap(seq.generalizedKmer);
	SEQAN_ASSERT_EQ(length(seq.generalizedKmer),length(back.generalizedKmer));
	SEQAN_ASSERT_EQ(length(seq.generalizedKmer),9);
	/*
		ACG --> WCG,MCG,RCG
				AYG,AMG,ASG
				ACK,ACR,ACS
	*/

	
	std::map<String<Iupac>,unsigned int > testMap;
	testMap["WCG"]=2;
	testMap["MCG"]=2;
	testMap["RCG"]=2;
	testMap["AYG"]=2;
	testMap["AMG"]=2;
	testMap["ASG"]=2;
	testMap["ACK"]=2;
	testMap["ACR"]=2;
	testMap["ACS"]=3;
	std::map<String<Iupac>,unsigned int >::iterator MapIterator;
	std::map<String<Iupac>,unsigned int >::iterator MapIterator2;
	
	MapIterator = seq.generalizedKmer.begin();
	MapIterator2=testMap.begin();
	SEQAN_ASSERT_EQ(length(seq.generalizedKmer), length(testMap));
	for(;MapIterator != seq.generalizedKmer.end();++MapIterator,++MapIterator2){
		std::cout<<(*MapIterator).first<<" "<<(*MapIterator2).first<<std::endl;
		SEQAN_ASSERT_EQ((*MapIterator).first,(*MapIterator2).first);
		SEQAN_ASSERT_EQ((*MapIterator).second,(*MapIterator2).second);

	}
	std::map<String<Iupac>,unsigned int> generalizedKmerTemp;
	std::map<String<Iupac>,unsigned int> generalizedKmerBackgroundTemp;
	generalizedKmerTemp =seq.generalizedKmer;
	generalizedKmerBackgroundTemp=back.generalizedKmer;
	seq.generalizedKmer.clear();
	back.generalizedKmer.clear();
	//Kmer zum r√ºbergehen, Temp zum neubef√ºllen
	GeneralizeKmer("ACS",generalizedKmerTemp,generalizedKmerBackgroundTemp, IMaps, seq, back);
	/*
		ACS -->ACB, ACV, AMS, AYS, ASS, RCS, WCS, MCS
	*/
	SEQAN_ASSERT_EQ(length(generalizedKmerTemp),length(generalizedKmerBackgroundTemp));
	SEQAN_ASSERT_EQ(length(seq.generalizedKmer),8);// ACB,ACV,AMS,AYS,ASS,RCS,WCS,MCS
	SEQAN_ASSERT_EQ(length(generalizedKmerTemp),9);// WCG,MCG,RCG,AYG,AMG,ASG,ACK,ACR,ACS
	//PrintMap(generalizedKmerTemp);

	generalizedKmerTemp=seq.generalizedKmer;
	generalizedKmerBackgroundTemp=back.generalizedKmer;
	seq.generalizedKmer.clear();
	back.generalizedKmer.clear();
	//umgekehrt
	GeneralizeKmer("ACB",generalizedKmerTemp,generalizedKmerBackgroundTemp,IMaps, seq, back);
	//PrintMap(seq.generalizedKmer);

	
	clear(seq.seqs);
	clear(back.seqs);
	clear(seq.SArray);
	clear(back.SArray);
	clear(seq.generalizedKmer);
	clear(back.generalizedKmer);
	clear(generalizedKmerTemp);
	clear(generalizedKmerBackgroundTemp);
	

}

SEQAN_DEFINE_TEST(test_ChIPSeq_initGeneralize){

	using namespace seqan;
	Seq seq;
	Seq back;

	seq.c=1;
	back.c=1;
	seq.seed=100;
	
	
	StringSet<String<Dna5> > test;
	appendValue(test,"ACGT");
	appendValue(test,"AGGT");
	appendValue(test,"ACCT");
	appendValue(test,"");

	StringSet<String<Dna5> > test2;
	appendValue(test2,"CCGT");
	appendValue(test2,"TTTT");
	appendValue(test2,"AACT");

	seq.seqs=test;
	back.seqs=test2;

	//typedef std::map<String<Dna5>,std::vector<int> > DnaCounterMap;
	//DnaCounterMap seqCount;
	//DnaCounterMap backCount;

	seq.SeqsNumber=length(seq.seqs);
	back.SeqsNumber=length(back.seqs);

	initST(seq);
	initST(back);
	initExactKmer(seq,back,3,3);
	//PrintMap(seq.seqCounter,seq.SeqsNumber);
	//std::cout<<std::endl;
	//PrintMap(back.seqCounter,back.SeqsNumber);
	/*std::map<String<Iupac>,double > generalizedSortedPValueReversed;
	typedef std::multimap<double,String<Dna5> > pValueMap;
	pValueMap SortedPValue;*/
	logFactorial(seq.SeqsNumber+back.SeqsNumber);//save all relevant factorial numbers
	FisherExactTest(seq,back);//computes the pValue of each Motif due to the counter
	//std::cout<<std::endl;
	PrintMap(seq.SortedPValue);
	std::multimap<double,String<Dna5> >::iterator SortedPValueIterator;
	SortedPValueIterator = seq.SortedPValue.begin();

	
	SEQAN_ASSERT_EQ(double(int(((*SortedPValueIterator).first+0.005)*100))/100,0.57);
	SEQAN_ASSERT_EQ((*SortedPValueIterator).second,"ACC");
	++SortedPValueIterator;
	SEQAN_ASSERT_EQ(double(int(((*SortedPValueIterator).first+0.005)*100))/100,0.57);
	SEQAN_ASSERT_EQ((*SortedPValueIterator).second,"ACG");
	++SortedPValueIterator;
	SEQAN_ASSERT_EQ(double(int(((*SortedPValueIterator).first+0.005)*100))/100,0.57);
	SEQAN_ASSERT_EQ((*SortedPValueIterator).second,"AGG");
	++SortedPValueIterator;
	SEQAN_ASSERT_EQ(double(int(((*SortedPValueIterator).first+0.005)*100))/100,0.57);
	SEQAN_ASSERT_EQ((*SortedPValueIterator).second,"CCT");
	++SortedPValueIterator;
	SEQAN_ASSERT_EQ(double(int(((*SortedPValueIterator).first+0.005)*100))/100,0.57);
	SEQAN_ASSERT_EQ((*SortedPValueIterator).second,"GGT");
	++SortedPValueIterator;
	SEQAN_ASSERT_EQ(double(int(((*SortedPValueIterator).first+0.005)*100))/100,0.86);
	SEQAN_ASSERT_EQ((*SortedPValueIterator).second,"CGT");
	++SortedPValueIterator;
	/***
		Teil aus InitGeneralization
	***/
	IupacMaps IMaps;
	MapIupac(IMaps);
	/*std::multimap<double,String<Dna5> >::iterator MapIterator;
	unsigned int i=0;
	unsigned int limit;
	if(seq.SortedPValue.size()>seq.seed)	limit=seq.seed;//seed meist = 100
	else if(seq.SortedPValue.size()==0) return;
	else	limit = seq.SortedPValue.size();
	for(MapIterator=seq.SortedPValue.begin();i<limit;++MapIterator,++i){//iterate over Top100
		GeneralizeKmer((*MapIterator).second,IMaps,seq,back);
	}



	SEQAN_ASSERT_EQ(seq.generalizedKmer.size(),back.generalizedKmer.size());
	//PrintMap(seq.generalizedKmer);
	SEQAN_ASSERT_EQ(seq.generalizedKmer.size(),50);
	SEQAN_ASSERT_EQ(seq.generalizedKmer["ACS"],2);
	SEQAN_ASSERT_EQ(seq.generalizedKmer["SCT"],1);
	SEQAN_ASSERT_EQ(seq.generalizedKmer["SGT"],2);*/
	//std::map<String<Iupac>,unsigned int> generalizedKmer;//unsigned int = estimated counter
	//std::map<String<Iupac>,unsigned int> generalizedKmerBackground;
	//std::multimap<double,String<Iupac> > generalizedSortedPValue;
	
	InitGeneralization(IMaps,seq,back);
	std::cout<<seq.generalizedSortedPValue.size()<<std::endl;
	PrintMap(seq.generalizedSortedPValue);


	clear(seq.seqs);
	clear(back.seqs);
	clear(seq.SArray);
	clear(back.SArray);
	clear(seq.generalizedKmer);
	clear(back.generalizedKmer);
	seq.generalizedSortedPValue.clear();

}

SEQAN_DEFINE_TEST(test_ChIPSeq_loopOverKmer){
	using namespace seqan;
	Seq seq;
	Seq back;

}


SEQAN_DEFINE_TEST(test_ChIPSeq_CounterGeneralize){

	using namespace seqan;
	Seq seq;
	Seq back;
	seq.seed=100;
	seq.c=1;
	back.c=1;
	
	StringSet<String<Dna5> > test;
	appendValue(test,"ACGACGTAGG");
	appendValue(test,"AGGT");
	appendValue(test,"AAGTGGG");
	appendValue(test,"AC");
	appendValue(test,"");

	StringSet<String<Dna5> > test2;
	appendValue(test2,"AAACCCGGG");
	appendValue(test2,"ATATATATAT");
	appendValue(test2,"CCTAAT");
	appendValue(test2,"GT");
	appendValue(test2,"");


	seq.seqs=test;
	back.seqs=test2;
	seq.SeqsNumber=length(seq.seqs);
	back.SeqsNumber=length(back.seqs);
	initST(seq);
	initST(back);

	//std::map<unsigned int,char> IupacMap;
	//std::map<char,unsigned int> IupacMapReversed;
	//std::map<char,String<Iupac> > IupacMapReplace; //stores the replacement-chars
	//std::map<char,String<Dna5> > IupacMapReplaceReversed;
	//MapIupac(IupacMap, IupacMapReversed, IupacMapReplace,IupacMapReplaceReversed);//IupacMap for generalization
	IupacMaps IMaps;
	MapIupac(IMaps);
	std::map<String<Iupac>,unsigned int > seqCounter;
	std::map<String<Iupac>,unsigned int > backCounter;
	Finder<Index<StringSet<String<Dna5> > > > finder(seq.SArray);
	Finder<Index<StringSet<String<Dna5> > > > finderB(back.SArray);//finder background
	String<Iupac> Kmer = "ASG";
	
	CountKmer(seqCounter,finder,Kmer,seq,IMaps);
	CountKmer(backCounter,finderB,Kmer,back,IMaps);
	Kmer="SGT";
	CountKmer(seqCounter,finder,Kmer,seq,IMaps);
	CountKmer(backCounter,finderB,Kmer,back,IMaps);
	Kmer="RGG";
	CountKmer(seqCounter,finder,Kmer,seq,IMaps);
	CountKmer(backCounter,finderB,Kmer,back,IMaps);
	Kmer="GGK";
	CountKmer(seqCounter,finder,Kmer,seq,IMaps);
	CountKmer(backCounter,finderB,Kmer,back,IMaps);
	Kmer="NGT";
	CountKmer(seqCounter,finder,Kmer,seq,IMaps);
	CountKmer(backCounter,finderB,Kmer,back,IMaps);
	Kmer="AAA";
	//std::cout<<Kmer;
	CountKmer(seqCounter,finder,Kmer,seq,IMaps);
	CountKmer(backCounter,finderB,Kmer,back,IMaps);
	Kmer="ACGACGTAGR";
	//std::cout<<Kmer;
	CountKmer(seqCounter,finder,Kmer,seq,IMaps);
	CountKmer(backCounter,finderB,Kmer,back,IMaps);
	Kmer="";
	CountKmer(seqCounter,finder,Kmer,seq,IMaps);
	CountKmer(backCounter,finderB,Kmer,back,IMaps);
	//std::cout<<std::endl;
	//std::cout<<std::endl;
	//std::cout<<Kmer;
	//PrintMap(seqCounter,seq.SeqsNumber);
	//std::cout<<std::endl;
	//PrintMap(backCounter,back.SeqsNumber);


	
		
	
	SEQAN_ASSERT_EQ(seqCounter[""],0);//should not appear in real data--> at least one wildcard
	SEQAN_ASSERT_EQ(seqCounter["AAA"],0);//should not appear in real data--> at least one wildcard
	SEQAN_ASSERT_EQ(seqCounter["ASG"],2);
	SEQAN_ASSERT_EQ(seqCounter["SGT"],2);
	SEQAN_ASSERT_EQ(seqCounter["RGG"],3);
	SEQAN_ASSERT_EQ(seqCounter["GGK"],2);
	SEQAN_ASSERT_EQ(seqCounter["NGT"],3);
	SEQAN_ASSERT_EQ(seqCounter["ACGACGTAGR"],1);

	
	SEQAN_ASSERT_EQ(backCounter[""],0);
	SEQAN_ASSERT_EQ(backCounter["AAA"],0);
	SEQAN_ASSERT_EQ(backCounter["ASG"],0);
	SEQAN_ASSERT_EQ(backCounter["SGT"],0);
	SEQAN_ASSERT_EQ(backCounter["RGG"],1);
	SEQAN_ASSERT_EQ(backCounter["GGK"],1);
	SEQAN_ASSERT_EQ(backCounter["NGT"],0);
	SEQAN_ASSERT_EQ(backCounter["ACGACGTAGR"],0);

	clear(seq.seqs);
	clear(back.seqs);
	clear(seq.SArray);
	clear(back.SArray);
	clear(seq.generalizedKmer);
	clear(back.generalizedKmer);
	seq.generalizedSortedPValue.clear();

}

SEQAN_DEFINE_TEST(test_ChIPSeq_estimateCounter){

	using namespace seqan;
	Seq seq;
	seq.seed=100;
	unsigned x = 10 + 10 - (10*10/10);
	SEQAN_ASSERT_EQ(x,10);
	x = 9 + 8 - (9*8/10);
	SEQAN_ASSERT_EQ(x,10);

	seq.SeqsNumber=11;
	seq.seqCounter["AAA"]=10;
	seq.seqCounter["TAA"]=3;
	seq.seqCounter["TGA"]=4;
	seq.seqCounter["ATA"]=9;
	seq.seqCounter["CAA"]=11;

	String<Iupac> temp = "AAA";
	String<Iupac> temp2= "TAA";
	estimateCounter(seq,temp,temp2,x);
	SEQAN_ASSERT_EQ(x,11);

	temp = "AAA";
	temp2= "CAA";
	estimateCounter(seq,temp,temp2,x);
	SEQAN_ASSERT_EQ(x,11);

	temp = "AAA";
	temp2= "ACA";
	estimateCounter(seq,temp,temp2,x);
	SEQAN_ASSERT_EQ(x,10);

	temp = "TAA";
	temp2= "TGA";
	estimateCounter(seq,temp,temp2,x);
	SEQAN_ASSERT_EQ(x,6);

	std::map<String<Iupac>,unsigned int> generalizedKmer;

	generalizedKmer["WAA"]=11;
	generalizedKmer["KAA"]=5;
	generalizedKmer["SAY"]=1;
	generalizedKmer["SAA"]=9;

	temp = "WAA";
	 temp2= "CAA";
	estimateCounter(seq,generalizedKmer,temp,temp2,x);
	SEQAN_ASSERT_EQ(x,11);

	temp = "KAA";
	temp2= "TTT";
	estimateCounter(seq,generalizedKmer,temp,temp2,x);
	SEQAN_ASSERT_EQ(x,5);

	temp = "SAA";
	temp2= "SAY";
	estimateCounter(seq,generalizedKmer,temp,temp2,x);
	SEQAN_ASSERT_EQ(x,10);

	clear(generalizedKmer);
	clear(seq.seqCounter);

}




SEQAN_BEGIN_TESTSUITE(dreme_test)
{
    // Call tests.
    SEQAN_CALL_TEST(test_ChIPSeq_fasta_print);
    SEQAN_CALL_TEST(test_ChIPSeq_suffixArray);
    SEQAN_CALL_TEST(test_ChIPSeq_initExactKmer);
    SEQAN_CALL_TEST(test_ChIPSeq_fact);
    SEQAN_CALL_TEST(test_ChIPSeq_modifyFET);
    SEQAN_CALL_TEST(test_ChIPSeq_generalize);
    SEQAN_CALL_TEST(test_ChIPSeq_initGeneralize);
    SEQAN_CALL_TEST(test_ChIPSeq_CounterGeneralize);
	SEQAN_CALL_TEST(test_ChIPSeq_estimateCounter);
}
SEQAN_END_TESTSUITE
