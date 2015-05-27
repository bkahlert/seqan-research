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


	//typedef Index<StringSet<String<Dna5> > > TMyIndex;//Dna5Q
	//Iterator<TMyIndex, BottomUp<> >::Type myIterator(seq.SArray);
	//for(;!atEnd(myIterator);++myIterator){
	//	std::cout<<representative(myIterator)<<std::endl;
	//	
	//}
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
	//PrintMap(seq.generalizedKmer);
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
	//	std::cout<<(*MapIterator).first<<" "<<(*MapIterator2).first<<std::endl;
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
	//PrintMap(seq.SortedPValue);
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
	std::multimap<double,String<Dna5> >::iterator MapIterator;
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
	SEQAN_ASSERT_EQ(seq.generalizedKmer["SGT"],2);
	//std::map<String<Iupac>,unsigned int> generalizedKmer;//unsigned int = estimated counter
	//std::map<String<Iupac>,unsigned int> generalizedKmerBackground;
	//std::multimap<double,String<Iupac> > generalizedSortedPValue;
	/*
	InitGeneralization(IMaps,seq,back);
	std::cout<<seq.generalizedSortedPValue.size()<<std::endl;
	PrintMap(seq.generalizedSortedPValue);
	*/

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

SEQAN_DEFINE_TEST(test_ChIPSeq_BuildFrequency){

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
	appendValue(test,"ACGAG");
	appendValue(test,"");

	StringSet<String<Dna5> > test2;
	appendValue(test2,"AAACCCGGG");
	appendValue(test2,"ATATATATAT");
	appendValue(test2,"CCTAAT");
	appendValue(test2,"GTGTGT");
	appendValue(test2,"");


	seq.seqs=test;
	back.seqs=test2;
	seq.SeqsNumber=length(seq.seqs);
	back.SeqsNumber=length(back.seqs);
	initST(seq);
	initST(back);
	priorFreq(seq);
	priorFreq(back);
	IupacMaps IMaps;
	MapIupac(IMaps);
	resize(seq.intervals, seq.SeqsNumber);
	resize(seq.intervalTrees, seq.SeqsNumber);

	
	Finder<Index<StringSet<String<Dna5> > > > finder(seq.SArray);
	String<Iupac>Kmer ="ASG";
	BuildFrequencyMatrix(finder,Kmer,seq,IMaps);//ASG gets masked

	SEQAN_ASSERT_EQ(double(int((seq.freqMatrix[0]['A']+0.005)*100))/100,0.86);
	SEQAN_ASSERT_EQ(seq.freqMatrix[0]['G'],0);
	SEQAN_ASSERT_EQ(seq.freqMatrix[1]['A'],0);
	SEQAN_ASSERT_EQ(double(int((seq.freqMatrix[1]['G']+0.005)*100))/100,0.5);//in der hälfte der sequenzen kommt ACG und AGG vor, ACG und AGG können in der selben sequenz sein
	SEQAN_ASSERT_EQ(double(int((seq.freqMatrix[1]['C']+0.005)*100))/100,0.5);
	SEQAN_ASSERT_EQ(seq.freqMatrix[2]['G'],1);
	
	for(unsigned s=0;s<seq.SeqsNumber;++s){
				createIntervalTree(seq.intervalTrees[s], seq.intervals[s]);
			}
	
	
	appendValue(test,"ACGACGTAGG");
	appendValue(test,"AGGT");
	appendValue(test,"AAGTGGG");
	appendValue(test,"ACGAG");
	++seq.c;
	initExactKmer(seq,back,3,3);
	//PrintMap(seq.seqCounter,seq.SeqsNumber);
	SEQAN_ASSERT_EQ(seq.seqCounter["AAG"],1);
	SEQAN_ASSERT_EQ(seq.seqCounter["ACG"],0);
	SEQAN_ASSERT_EQ(seq.seqCounter["AGT"],1);
	SEQAN_ASSERT_EQ(seq.seqCounter["GTG"],1);
	SEQAN_ASSERT_EQ(seq.seqCounter["GGG"],1);
}

SEQAN_DEFINE_TEST(test_PWMClustering){

	using namespace seqan;

	
	std::vector<std::vector<String<double>> > compare(5, std::vector<String<double>>(5));
	String<double> M;
	appendValue(M,0);
	String<double> minDifference;
	resize(minDifference,3);
	minDifference[0]=0;
	/*****
			Matrix:
			0
			2	0
			3	3	0
			4	1	7	0
			6	4	6	9	0


	*****/

	compare[0][0]=M;
	M[0]=2;
	compare[1][0]=M;
	M[0]=0;
	compare[1][1]=M;
	M[0]=3;
	compare[2][0]=M;
	M[0]=3;
	compare[2][1]=M;
	M[0]=0;
	compare[2][2]=M;
	M[0]=4;
	compare[3][0]=M;
	M[0]=1;
	compare[3][1]=M;
	M[0]=7;
	compare[3][2]=M;
	M[0]=0;
	compare[3][3]=M;
	M[0]=6;
	compare[4][0]=M;
	M[0]=4;
	compare[4][1]=M;
	M[0]=6;
	compare[4][2]=M;
	M[0]=9;
	compare[4][3]=M;
	M[0]=0;
	compare[4][4]=M;
	

	//minDifference[0]=compare[3][1][0];//3=zeile 1=spalte
	//minDifference[1]=1;
	//minDifference[2]=3;
	unsigned j;
	unsigned n=5;
	minDifferenceInMatrix(n,minDifference,compare);
	
	UpdateDistantMatrix(n,int(minDifference[1]),int(minDifference[2]),compare,CompleteLinkage());


	
	SEQAN_ASSERT_EQ(compare[0][0][0],0);
	SEQAN_ASSERT_EQ(compare[1][0][0],4);
	SEQAN_ASSERT_EQ(compare[1][1][0],0);
	SEQAN_ASSERT_EQ(compare[2][0][0],3);
	SEQAN_ASSERT_EQ(compare[2][1][0],7);
	SEQAN_ASSERT_EQ(compare[2][2][0],0);
	SEQAN_ASSERT_EQ(compare[3][0][0],6);
	SEQAN_ASSERT_EQ(compare[3][1][0],9);
	SEQAN_ASSERT_EQ(compare[3][2][0],6);
	SEQAN_ASSERT_EQ(compare[3][3][0],0);
	/*****
			New Matrix:

			0
			4	0
			3	7	0
			6	9	6	0
	*****/

	//minDifference[0]=compare[2][0][0];//3=zeile 1=spalte
	//minDifference[1]=0;
	//minDifference[2]=2;
	--n;
	minDifferenceInMatrix(n,minDifference,compare);
	UpdateDistantMatrix(n,int(minDifference[1]),int(minDifference[2]),compare,CompleteLinkage());

	/*****
			New Matrix:

			0
			7	0
			6	9	0
			
	*****/
	SEQAN_ASSERT_EQ(compare[0][0][0],0);
	SEQAN_ASSERT_EQ(compare[1][0][0],7);
	SEQAN_ASSERT_EQ(compare[1][1][0],0);
	SEQAN_ASSERT_EQ(compare[2][0][0],6);
	SEQAN_ASSERT_EQ(compare[2][1][0],9);
	SEQAN_ASSERT_EQ(compare[2][2][0],0);




}

SEQAN_DEFINE_TEST(test_ComparePWM){

	using namespace seqan;

	std::map<Iupac,double>  freqMatrix1;
	std::map<Iupac,double>  freqMatrix2;
	double test=0;

	freqMatrix1['A']=0.25;
	freqMatrix1['C']=0.25;
	freqMatrix1['G']=0.25;
	freqMatrix1['T']=0.25;

	freqMatrix2['A']=0.97;
	freqMatrix2['C']=0.01;
	freqMatrix2['G']=0.01;
	freqMatrix2['T']=0.01;

	test=double(int((ComparePWM(freqMatrix1,freqMatrix2,Entropy())+0.005)*100))/100;
	SEQAN_ASSERT_EQ(test,1.65);
	


	freqMatrix1['A']=0.03;
	freqMatrix1['C']=0.75;
	freqMatrix1['G']=0.03;
	freqMatrix1['T']=0.03;

	freqMatrix2['A']=0.3;
	freqMatrix2['C']=0.01;
	freqMatrix2['G']=0.33;
	freqMatrix2['T']=0.36;

	test=double(int((ComparePWM(freqMatrix1,freqMatrix2,Entropy())+0.005)*100))/100;
	SEQAN_ASSERT_EQ(test,2.68);

	freqMatrix1['A']=0.3;
	freqMatrix1['C']=0.3;
	freqMatrix1['G']=0.3;
	freqMatrix1['T']=0.1;

	freqMatrix2['A']=0.3;
	freqMatrix2['C']=0.3;
	freqMatrix2['G']=0.3;
	freqMatrix2['T']=0.1;

	test=double(int((ComparePWM(freqMatrix1,freqMatrix2,Entropy())+0.005)*100))/100;
	SEQAN_ASSERT_EQ(test,0);

	freqMatrix1['A']=0.4;
	freqMatrix1['C']=0.15;
	freqMatrix1['G']=0.4;
	freqMatrix1['T']=0.05;

	freqMatrix2['A']=0.5;
	freqMatrix2['C']=0.1;
	freqMatrix2['G']=0.3;
	freqMatrix2['T']=0.1;

	test=double(int((ComparePWM(freqMatrix1,freqMatrix2,Entropy())+0.005)*100))/100;
	SEQAN_ASSERT_EQ(test,0.05);

}


SEQAN_DEFINE_TEST(test_AlignPWMs){

	using namespace seqan;

	std::map<unsigned int,std::map<Iupac,double> >  freqMatrix1;
	std::map<unsigned int,std::map<Iupac,double> >  freqMatrix2;
	String<double> Mmax;

	freqMatrix1[0]['A']=0.66;
	freqMatrix1[0]['C']=0.3;
	freqMatrix1[0]['G']=0.02;
	freqMatrix1[0]['T']=0.02;

	freqMatrix1[1]['A']=0.65;
	freqMatrix1[1]['C']=0.025;
	freqMatrix1[1]['G']=0.3;
	freqMatrix1[1]['T']=0.025;

	freqMatrix1[2]['A']=0.3;
	freqMatrix1[2]['C']=0.3;
	freqMatrix1[2]['G']=0.3;
	freqMatrix1[2]['T']=0.1;

	freqMatrix1[3]['A']=0.45;
	freqMatrix1[3]['C']=0.05;
	freqMatrix1[3]['G']=0.05;
	freqMatrix1[3]['T']=0.45;

	freqMatrix1[4]['A']=0.03;
	freqMatrix1[4]['C']=0.03;
	freqMatrix1[4]['G']=0.04;
	freqMatrix1[4]['T']=0.9;
	/*********************/
	freqMatrix2[0]['A']=0.3;
	freqMatrix2[0]['C']=0.05;
	freqMatrix2[0]['G']=0.6;
	freqMatrix2[0]['T']=0.05;

	freqMatrix2[1]['A']=0.4;
	freqMatrix2[1]['C']=0.1;
	freqMatrix2[1]['G']=0.1;
	freqMatrix2[1]['T']=0.4;

	freqMatrix2[2]['A']=0.33;
	freqMatrix2[2]['C']=0.33;
	freqMatrix2[2]['G']=0.33;
	freqMatrix2[2]['T']=0.01;

	

	Mmax=AlignPWMs(freqMatrix1,freqMatrix2);
	SEQAN_ASSERT_EQ(Mmax[1],4);
	SEQAN_ASSERT_EQ(Mmax[2],3);//das 3-mer sitzt mitten im 5-mer --> Überlappung der Stellen 2,3,4 im 5-mer
	
	clear(freqMatrix1);
	clear(freqMatrix2);

	freqMatrix1[0]['A']=0.66;
	freqMatrix1[0]['C']=0.3;
	freqMatrix1[0]['G']=0.02;
	freqMatrix1[0]['T']=0.02;

	freqMatrix1[1]['A']=0.65;
	freqMatrix1[1]['C']=0.025;
	freqMatrix1[1]['G']=0.3;
	freqMatrix1[1]['T']=0.025;

	freqMatrix1[2]['A']=0.3;
	freqMatrix1[2]['C']=0.3;
	freqMatrix1[2]['G']=0.3;
	freqMatrix1[2]['T']=0.1;

	

	
	/*********************/
	freqMatrix2[0]['A']=0.45;
	freqMatrix2[0]['C']=0.45;
	freqMatrix2[0]['G']=0.05;
	freqMatrix2[0]['T']=0.05;

	freqMatrix2[1]['A']=0.9;
	freqMatrix2[1]['C']=0.033;
	freqMatrix2[1]['G']=0.033;
	freqMatrix2[1]['T']=0.034;

	freqMatrix2[2]['A']=0.33;
	freqMatrix2[2]['C']=0.33;
	freqMatrix2[2]['G']=0.33;
	freqMatrix2[2]['T']=0.01;

	freqMatrix2[3]['A']=0.05;
	freqMatrix2[3]['C']=0.05;
	freqMatrix2[3]['G']=0.05;
	freqMatrix2[3]['T']=0.85;

	freqMatrix2[4]['A']=0.45;
	freqMatrix2[4]['C']=0.05;
	freqMatrix2[4]['G']=0.05;
	freqMatrix2[4]['T']=0.45;
	
	Mmax=AlignPWMs(freqMatrix1,freqMatrix2);
	
	SEQAN_ASSERT_EQ(Mmax[1],3);
	SEQAN_ASSERT_EQ(Mmax[2],3);//das 3-mer sitzt mitten im 5-mer --> Überlappung der Stellen 2,3,4 im 5-mer

}

SEQAN_DEFINE_TEST(test_computeDr){

	using namespace seqan;

	String<int> temp;
	Cluster cluster;
	appendValue(temp,1,Exact());
	appendValue(temp,2,Exact());
	appendValue(temp,4,Exact());
	cluster.content=temp;
	clear(temp);

	std::vector<std::vector<String<double> > > compare(5, std::vector<String<double> >(5));
	compare[1][0]=1.2;
	compare[2][0]=0.9;
	compare[3][0]=1.34;
	compare[4][0]=2.67;
	compare[2][1]=0.75;
	compare[3][1]=0.45;
	compare[4][1]=1.23;
	compare[3][2]=1.74;
	compare[4][2]=1.01;
	compare[4][3]=0.97;

	SEQAN_ASSERT_EQ(double(int((computeDr(cluster,compare)+0.005)*100))/100,2.99);
	clear(cluster.content);

	appendValue(temp,2,Exact());
	appendValue(temp,3,Exact());
	appendValue(temp,1,Exact());
	cluster.content=temp;
	
	SEQAN_ASSERT_EQ(double(int((computeDr(cluster,compare)+0.005)*100))/100,2.94);
	clear(cluster.content);
	clear(temp);

	clear(compare);

}

SEQAN_DEFINE_TEST(test_computeWk){

	using namespace seqan;

	String<Cluster> cluster;
	resize(cluster,7);
	int n=0;
	String<int> temp;

	cluster[n].left=1;
	cluster[n].right=2;
	appendValue(temp,1);
	appendValue(temp,2);
	cluster[n].content=temp;
	cluster[n].Dr=2.99;
	clear(temp);
	
	SEQAN_ASSERT_EQ(double(int((computeWk(n,cluster)+0.0005)*1000))/1000,1.495);
	SEQAN_ASSERT_EQ(double(int((cluster[n].Wk+0.0005)*1000))/1000,1.495);

	++n;
	cluster[n].left=-1;
	cluster[n].right=3;
	appendValue(temp,1);
	appendValue(temp,2);
	appendValue(temp,3);
	cluster[n].content=temp;
	cluster[n].Dr=3.33;
	SEQAN_ASSERT_EQ(double(int((computeWk(n,cluster)+0.005)*100))/100,1.11);
	SEQAN_ASSERT_EQ(double(int((cluster[n].Wk+0.005)*100))/100,1.11);
	clear(temp);

	++n;
	cluster[n].left=4;
	cluster[n].right=5;
	appendValue(temp,5);
	appendValue(temp,4);
	cluster[n].content=temp;
	cluster[n].Dr=2.5;
	SEQAN_ASSERT_EQ(double(int((computeWk(n,cluster)+0.005)*100))/100,2.36);
	SEQAN_ASSERT_EQ(double(int((cluster[n].Wk+0.005)*100))/100,2.36);
	clear(temp);

	++n;
	cluster[n].left=-2;
	cluster[n].right=-3;
	appendValue(temp,1);
	appendValue(temp,2);
	appendValue(temp,3);
	appendValue(temp,5);
	appendValue(temp,4);
	cluster[n].content=temp;
	cluster[n].Dr=4.5;
	SEQAN_ASSERT_EQ(double(int((computeWk(n,cluster)+0.05)*10))/10,0.9);
	SEQAN_ASSERT_EQ(double(int((cluster[n].Wk+0.05)*10))/10,0.9);
	clear(temp);
	clear(cluster);
}

SEQAN_DEFINE_TEST(test_compute_l_quer){

	using namespace seqan;

	int B=3;
	unsigned allPWMsLength = 4;

	String<Cluster> cluster;
	String<double> l_quer_for_k;
	String<String<double> > allWk;

	resize(allWk,B);
	resize(cluster,allPWMsLength);
	resize(l_quer_for_k,allPWMsLength);
	l_quer_for_k[0]=0;

	unsigned i=0;
	cluster[0].Wk=0.5;
	cluster[1].Wk=0.6;
	cluster[2].Wk=0.44;	
	compute_l_quer(cluster,l_quer_for_k,allWk,i,B,allPWMsLength);
	clear(cluster);
	resize(cluster,allPWMsLength);

	++i;
	cluster[0].Wk=0.9;
	cluster[1].Wk=0.2;
	cluster[2].Wk=0.11;	
	compute_l_quer(cluster,l_quer_for_k,allWk,i,B,allPWMsLength);
	clear(cluster);
	resize(cluster,allPWMsLength);

	++i;
	cluster[0].Wk=0.13;
	cluster[1].Wk=0.16;
	cluster[2].Wk=0.414;
	compute_l_quer(cluster,l_quer_for_k,allWk,i,B,allPWMsLength);
	clear(cluster);
	
	
	SEQAN_ASSERT_EQ(double(int((l_quer_for_k[1]-0.0005)*1000))/1000,-1.303);
	SEQAN_ASSERT_EQ(double(int((l_quer_for_k[2]-0.0005)*1000))/1000,-1.318);
	SEQAN_ASSERT_EQ(double(int((l_quer_for_k[3]-0.0005)*1000))/1000,-0.946);
	SEQAN_ASSERT_EQ(double(int((l_quer_for_k[0]-0.0005)*1000))/1000,-3.567);
	
}

SEQAN_DEFINE_TEST(test_computeGap){

	using namespace seqan;

	unsigned allPWMsLength=4;
	String<double> Gap;
	String<double> l_quer_for_k;
	String<Cluster> observedCluster;
	resize(observedCluster,allPWMsLength);

	resize(Gap,allPWMsLength);
	resize(l_quer_for_k,allPWMsLength);

	l_quer_for_k[1]=-1.303;
	l_quer_for_k[2]=-1.318;
	l_quer_for_k[3]=-0.946;

	observedCluster[0].Wk=0.976;
	observedCluster[1].Wk=0.364;
	observedCluster[2].Wk=0.782;

	computeGap(l_quer_for_k,observedCluster,Gap,allPWMsLength);

	SEQAN_ASSERT_EQ(double(int((Gap[1]-0.0005)*1000))/1000,-1.057);
	SEQAN_ASSERT_EQ(double(int((Gap[2]-0.0005)*1000))/1000,-0.307);
	SEQAN_ASSERT_EQ(double(int((Gap[3]-0.0005)*1000))/1000,-0.922);
}

SEQAN_DEFINE_TEST(test_compute_sdk_and_sk){

	using namespace seqan;

	String<double> l_quer_for_k;
	String<double> sk;
	String<String<double> > allWk;
	unsigned allPWMsLength=4;
	int B=3;
	resize(l_quer_for_k,allPWMsLength);	
	resize(sk,allPWMsLength);
	resize(allWk,B);

	l_quer_for_k[0]=-3.567;

	appendValue(allWk[0],0.5);
	appendValue(allWk[0],0.6);
	appendValue(allWk[0],0.44);
	
	appendValue(allWk[1],0.9);
	appendValue(allWk[1],0.2);
	appendValue(allWk[1],0.11);

	appendValue(allWk[2],0.13);
	appendValue(allWk[2],0.16);
	appendValue(allWk[2],0.414);


	compute_sdk_and_sk(sk,l_quer_for_k,allWk,allPWMsLength,B);

	SEQAN_ASSERT_EQ(double(int((sk[1]+0.0005)*1000))/1000,2.716);
	SEQAN_ASSERT_EQ(double(int((sk[2]+0.0005)*1000))/1000,2.682);
	SEQAN_ASSERT_EQ(double(int((sk[3]+0.0005)*1000))/1000,3.167);



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
	SEQAN_CALL_TEST(test_ChIPSeq_BuildFrequency);
	SEQAN_CALL_TEST(test_PWMClustering);
	SEQAN_CALL_TEST(test_ComparePWM);
	SEQAN_CALL_TEST(test_AlignPWMs);
	SEQAN_CALL_TEST(test_computeDr);
	SEQAN_CALL_TEST(test_computeWk);
	SEQAN_CALL_TEST(test_compute_l_quer);
	SEQAN_CALL_TEST(test_computeGap);
	SEQAN_CALL_TEST(test_compute_sdk_and_sk);
}
SEQAN_END_TESTSUITE
