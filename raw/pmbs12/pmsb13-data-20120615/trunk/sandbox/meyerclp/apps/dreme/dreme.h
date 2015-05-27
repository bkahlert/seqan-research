#ifndef SANDBOX_MEYERCLP_APPS_DREME_H_
#define SANDBOX_MEYERCLP_APPS_DREME_H_

#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/stream.h>
#include <seqan/find_motif.h>
#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/misc/misc_interval_tree.h>
#include <seqan/random.h>
#include <seqan/map.h>
//#include <seqan/seq_io.h>



using namespace seqan;


struct Seq
{
	
	unsigned c; 
	unsigned NumberOfMotifs;
	unsigned seed;
	unsigned SeqsNumber;
	double pValue;
	double treshold;

	StringSet<String<Dna5> > seqs;
	StringSet<CharString> ids;
	String<String<Iupac> > allMotifs;
	String<unsigned> results;
	String< std::map<unsigned int,std::map<Iupac,double> > > allPWMs;
	typedef IntervalAndCargo<unsigned, unsigned> TInterval;
	String<String<TInterval> > intervals;//String<TInterval> enthält alle Intervalle einer Sequenz. String<String<..>> enthält alle Sequenzen	
	String<IntervalTree<unsigned> > intervalTrees; // Tree für schnellere Suche erstellen --> String von Trees, da mehrere Sequenzen
	/**
		In BuildFrequencyMatrix werden die gefundenen Intervalle des Top-Motivs an Intervalls gehängt
		Am Ende der while wird dann ein IntervallTree für jede Sequenz erzeugt, sodass im nächsten Schritt direkt geschaut werden kann ob ein Motiv mit einem schon 
		gefundenen überlappt
	**/
	
	Index< StringSet<String<Dna5> > > SArray;
	
	std::map<String<Dna5>,unsigned int > seqCounter;//maps the Sequence-Kmere to a Counter for the Sequence
	std::map<String<Iupac>,unsigned int> generalizedKmer;
	std::map<String<Iupac>,double > SortedPValueReversed;
	std::multimap<double,String<Dna5> > SortedPValue;
	std::multimap<double,String<Iupac> > generalizedSortedPValue;
	std::map<unsigned int,std::map<Iupac,double> > freqMatrix;
	std::map<unsigned int,std::map<Iupac,double> > weightMatrix;
	std::map<unsigned int, double > seqLogoMatrix;
	std::map<unsigned int, double > columnEntropy;
	std::map<unsigned int,std::map<Iupac,double> > InformationContentMatrix;
	FrequencyDistribution<Dna5> frequencies;

	
	
	
};

struct IupacMaps
{

	std::map<unsigned int,char> IupacMap;
	std::map<char,unsigned int> IupacMapReversed;
	std::map<char,String<Iupac> > IupacMapReplace; //stores the replacement-chars
	std::map<char,String<Dna5> > IupacMapReplaceReversed;
	std::map<String<Dna5>, char > IupacMapInversed;
};


struct Cluster
{
	
	int left;
	int right;
	String<int>  content;
	double Dr;
	double Wk;


};


struct Euclid_;
typedef Tag<Euclid_> Euclid;

struct Pearson_;
typedef Tag<Pearson_> Pearson;

struct Entropy_;
typedef Tag<Entropy_> Entropy;

struct CompleteLinkage_;
typedef Tag<CompleteLinkage_> CompleteLinkage;

struct SingleLinkage_;
typedef Tag<SingleLinkage_> SingleLinkage;

struct AverageLinkage_;
typedef Tag<AverageLinkage_> AverageLinkage;

void readFastA(struct Seq &seq, CharString fname);
template <typename TStream>
void PrintFastA(TStream & stream, Seq &seq);
void initST(Seq &seq);
void PrintST(Seq &seq);
void computeSettings(Seq &sequences, unsigned &kmer_len, unsigned &kmer_len_end, bool &save, bool &clusterData);
void priorFreq(Seq &seq);
void initExactKmer(Seq &seq,Seq &back,unsigned int kmer_len,unsigned int kmer_len_end);
void CountKmer(Seq &seq, Finder<Index<StringSet<String<Dna5> > > > &finder, String<Dna5> &Kmer);
void CountKmer(std::map<String<Iupac>,unsigned int > &Dna5CounterMap, Finder<Index<StringSet<String<Iupac> > > > &finder, String<Iupac> &Kmer,Seq &seq,IupacMaps &IMap);
void PrintMap(std::map<String<Dna5>,unsigned int > &Dna5CounterMap,unsigned int SeqsNumber);
void PrintMap(std::map<String<Iupac>,unsigned int > &Dna5CounterMap,unsigned int SeqsNumber);
void PrintMap(std::multimap<double,String<Dna5> > &pValueMap);
void PrintMap(std::map<String<Iupac>,unsigned int> &generalizedKmer);
void PrintMap(Seq &seq,bool foreground);
void DebugMap(Seq &seq,Seq &back,std::map<String<Dna5>,unsigned int > &sequencesCounter,std::map<String<Dna5>,unsigned int > &backgroundCounter);
void DebugMultiMap(std::map<String<Dna5>,unsigned int > &sequencesCounter,std::multimap<double,String<Dna5> > &SortedPValue);
double* logFac;
void logFactorial(unsigned int len);
double calcFET(unsigned int a,unsigned int b,unsigned int c,unsigned int d);
void modifyFET(unsigned int a,unsigned int b,unsigned int c,unsigned int d, double &pValue);
void DebugFisherExact(unsigned int a,unsigned int b,unsigned int c,unsigned int d);
void FisherExactTest(Seq &seq, Seq &back);

void FisherExactTest(std::map<String<Iupac>,unsigned int > &SequenceCounter,std::map<String<Iupac>,unsigned int > &BackgroundCounter, Seq &seq, Seq &back);
double FisherExactTest(Seq &seq, Seq &back,std::multimap<double,String<Dna5> > &GeneralizedSortedPValueTemp);
void MapIupac(IupacMaps &IMaps);

void InitGeneralization(IupacMaps &IMaps,Seq &seq, Seq &back);
void loopOverKmer(Seq &seq,String<Iupac> &temp,String<Iupac> &Kmer,Iterator<String<Iupac> >::Type &tempIt,Finder<Index<StringSet<String<Dna5> > > > &finder,unsigned int &counter,std::vector<int> &CounterV,IupacMaps &IMap);
void FindKmer(Seq &seq,String<Iupac> &temp,Finder<Index<StringSet<String<Dna5> > > > &finder,unsigned int &counter,std::vector<int> &CounterV);
void GeneralizeKmer(String<Dna5> Kmer, IupacMaps &IMaps,Seq &seq, Seq &back);
void GeneralizeKmer(String<Iupac> Kmer,std::map<String<Iupac>,unsigned int> &generalizedKmerTemp,std::map<String<Iupac>,unsigned int> &generalizedKmerBackgroundTemp,IupacMaps &IMaps,Seq &seq, Seq &back);

void estimateCounter(Seq &seq,String<Iupac> temp,String<Iupac> temp2,unsigned int &counter);
void estimateCounter(Seq &seq,std::map<String<Iupac>,unsigned int> &generalizedKmer,String<Iupac> temp,String<Iupac> temp2,unsigned int &counter);
void exactGeneralizeCount(std::multimap<double,String<Iupac> > &SortedPValueG,std::map<String<Iupac>,unsigned int > &seqCounter,std::map<String<Iupac>,unsigned int > &backCounter,Finder<Index<StringSet<String<Dna5> > > > &finder,Finder<Index<StringSet<String<Dna5> > > > &finderB,Seq &seq, Seq &back,IupacMaps &IMap);
void FindTopKmer(Seq &seq,String<Iupac> &temp,Finder<Index<StringSet<String<Dna5> > > > &finder,unsigned int &counter,std::vector<int> &CounterV);
void FindTopKmer(Seq &seq,String<Dna5> &temp,Finder<Index<StringSet<String<Dna5> > > > &finder,unsigned int &counter,std::vector<int> &CounterV);
void loopOverTopKmer(	Seq &seq,String<Iupac> &temp,String<Iupac> &Kmer,Iterator<String<Iupac> >::Type &tempIt,Finder<Index<StringSet<String<Dna5> > > > &finder,unsigned int &counter,std::vector<int> &CounterV,IupacMaps &IMap);
void BuildFrequencyMatrix( Finder<Index<StringSet<String<Dna5> > > > &finder, String<Iupac> &Kmer,Seq &seq, IupacMaps &IMaps);
void BuildFrequencyMatrix( Finder<Index<StringSet<String<Dna5> > > > &finder,String<Dna5> &Kmer,Seq &seq, IupacMaps &IMaps);
void BuildWeightMatrix(Seq &seq);
void BuildInformationContentMatrix(Seq &seq);
double ComparePWM(std::map<Iupac,double>  &freqMatrix1,std::map<Iupac,double>  &freqMatrix2, Entropy const & tag);
double ComparePWM(std::map<Iupac,double>  &freqMatrix1,std::map<Iupac,double>  &freqMatrix2, Euclid const & tag);
double ComparePWM(std::map<Iupac,double>  &freqMatrix1,std::map<Iupac,double>  &freqMatrix2, Pearson const & tag);
void PWMClustering(String< std::map<unsigned int,std::map<Iupac,double> > > &allPWMs, String<Cluster> &cluster,String<int> &traceback);
void UpdateDistantMatrix(int n, int x, int y, std::vector<std::vector<String<double > > > &compare, CompleteLinkage const & tag);
void UpdateDistantMatrix(int n, int x, int y, std::vector<std::vector<String<double > > > &compare, SingleLinkage const & tag);
void UpdateDistantMatrix(int n, int x, int y, std::vector<std::vector<String<double > > > &compare,std::vector<unsigned> &weights, AverageLinkage const & tag);
void minDifferenceInMatrix(unsigned n,String<double> &minDifference,std::vector<std::vector<String<double > > > compare);
void replaceKmer(Seq &seq,unsigned int stringNumber, unsigned int begin, unsigned int end);
void saveData(Seq &seq,std::ofstream &PWM,unsigned c);
void computesDistantMatrix(String< std::map<unsigned int,std::map<Iupac,double> > > &allPWMs,std::vector<std::vector<String<double > > > &compare, unsigned allPWMsLength);
String<double> AlignPWMs(std::map<unsigned int,std::map<Iupac,double> > &freqMatrix1,std::map<unsigned int,std::map<Iupac,double> > &freqMatrix2);
double computeDr(Cluster &cluster,std::vector<std::vector<String<double> > > copyCompare);
double computeWk(int n,String<Cluster> &cluster);
void computeReferenceData(String< std::map<unsigned int,std::map<Iupac,double> > > &Reference,unsigned allPWMsLength,unsigned PWMLength );
unsigned computeGapStat(String< std::map<unsigned int,std::map<Iupac,double> > > &allPWMs, String<Cluster> &observedCluster);
void compute_l_quer(String<Cluster> &cluster, String<double> &l_quer_for_k,String<String<double> > &allWk,unsigned i,int B,	unsigned allPWMsLength);
void computeGap(String<double> &l_quer_for_k,String<Cluster> &observedCluster,String<double> &Gap,unsigned allPWMsLength);
void compute_sdk_and_sk(String<double> &sk,String<double> &l_quer_for_k,String<String<double> > &allWk,unsigned allPWMsLength,int B);
void clusterTraceback(String<String<Iupac> > &allMotifs,String<int> &traceback, unsigned j);


void readFastA(	struct Seq &seq, 
				CharString fname){

	
	std::fstream fasta(toCString(fname), std:: ios_base::in | std::ios_base::binary);
	if(!fasta.good()){
		std::cerr << "Could not open file: \"" << fname << "\"" << std::endl;
		std::exit(1);
	}
	typedef String<char,MMap<> > TMMapString;
	TMMapString mmapString;
	if(!open(mmapString, toCString(fname), OPEN_RDONLY))
		std::exit(1);
	RecordReader<TMMapString, DoublePass<Mapped> > reader(mmapString);


	AutoSeqStreamFormat formatTag;
	if(!checkStreamFormat(reader,formatTag)){
		std::cerr<<"Could not determine file format!"<<std::endl;
		std::exit(1);
	}
	std::cout<<"File format is "<<getAutoSeqStreamFormatName(formatTag)<<'\n';


	StringSet<CharString> ids;
	StringSet<String<Dna5> > seqs;//

	if(read2(ids,seqs,reader,formatTag) !=0){
				std::cerr<<"ERROR reading FASTA"<<std::endl;
				std::exit(1);
	}
	
	seq.ids = ids;
	seq.seqs = seqs;
}


template <typename TStream>
void PrintFastA(TStream & stream, 
				Seq &seq){

	SEQAN_ASSERT_EQ(length(seq.ids), length(seq.seqs));

	typedef Iterator<StringSet<CharString>, Rooted>::Type TIdIter;
	typedef Iterator<StringSet<String<Dna5> >, Standard>::Type TSeqIter;//Dna5Q
	TIdIter idIt =begin(seq.ids, Rooted());
	TSeqIter seqIt=begin(seq.seqs, Standard());
	for(;!atEnd(idIt);++idIt,++seqIt){

		
		stream <<*idIt<<'\t'<<*seqIt<<std::endl;

	}

	
}

void initST(Seq &seq){
	//creates a index based on an enhanced suffix array 
	indexText(seq.SArray) = seq.seqs;
}

void PrintST(Seq &seq){

	typedef Index<StringSet<String<Dna5> > > TMyIndex;
	Iterator<TMyIndex, BottomUp<> >::Type myIterator(seq.SArray);
	for(;!atEnd(myIterator);++myIterator){
		std::cout<<representative(myIterator)<<std::endl;
		
	}

}
/***
		Command-Line settings 
		ToDo: Error-Handling, Cluster-Settings
***/
void computeSettings(Seq &sequences, unsigned &kmer_len, unsigned &kmer_len_end, bool &save, bool &clusterData){
	char setting;
	std::cout<<"-------------DREME-Settings-------------"<<std::endl;
	std::cout<<"----------------------------------------"<<std::endl;
	std::cout<<"Standard-Settings are: mink(6), maxk(6), Seed(100), pValueTreshold(0.05), NumberOfMotifs(5)"<<std::endl;
	while(true){
		std::cout<<"Alter the Settings? (Y) (N)"<<std::endl;
		std::cin>>setting;
		if(setting =='Y' || setting=='N' ||setting =='y' ||setting =='n') break;
		std::cout<<"**** ERROR: MUST BE Y OR N ****"<<std::endl;
	}
	
	if(setting=='Y' || setting =='y'){
		std::cout<<"Enter in mink(should be >=3) "<<std::endl;
		std::cin>>kmer_len;
		std::cout<<"Enter in maxk(a value>8 could slow down the program) "<<std::endl;
		std::cin>>kmer_len_end;
		std::cout<<"Enter in Seed(100) "<<std::endl;
		std::cin>>sequences.seed;
		std::cout<<"Enter in pValueTreshold(0.05) "<<std::endl;
		std::cin>>sequences.treshold;
		std::cout<<"Enter in NumberOfMotifs(5) "<<std::endl;
		std::cin>>sequences.NumberOfMotifs;
		++sequences.NumberOfMotifs;


	}

	else{
		kmer_len=6;
		kmer_len_end=6;
		sequences.seed=100;
		sequences.treshold=0.05;
		sequences.NumberOfMotifs=6;
		
	}
	save=false;
	while(true){
		std::cout<<"Save the data (PWMs) to a file? (Y) (N)"<<std::endl;//path hardcoded--> ToDo
		std::cin>>setting;
		if(setting =='Y' || setting=='y'){
			save=true;
			break;
		}
		else if(setting =='N' ||setting =='n') break;
		std::cout<<"**** ERROR: MUST BE Y OR N ****"<<std::endl;
	}
	clusterData=false;
	while(true){
		std::cout<<"Should the data be clustered? (Y)(N)"<<std::endl;
		std::cin>>setting;
		if(setting =='Y' || setting=='y'){
			clusterData=true;
			break;
		}
		else if(setting =='N' ||setting =='n') break;
		std::cout<<"**** ERROR: MUST BE Y OR N ****"<<std::endl;

	}
}

/***
	Berechnet relative Hintergrund-Wahrscheinlichkeit von ACGT
***/
void priorFreq(Seq &seq){

	typedef Iterator<StringSet<String<Dna5> >, Standard>::Type TSeqIter;
	TSeqIter seqItbegin=begin(seq.seqs, Standard());
	TSeqIter seqItend=end(seq.seqs, Standard());
	absFreqOfLettersInSetOfSeqs(seq.frequencies,seqItbegin,seqItend);
	normalize(seq.frequencies);

}


//Initiiert Suche im Vorder- und Hintergrund
void initExactKmer( Seq &seq,
					Seq &back,
					unsigned int kmer_len,
					unsigned int kmer_len_end){

	Finder<Index<StringSet<String<Dna5> > > > finder(seq.SArray);
	Finder<Index<StringSet<String<Dna5> > > > finderB(back.SArray);
	typedef Index< StringSet<String<Dna5> > > TMyIndex;
	
	if(kmer_len<3) kmer_len=3;
	if(kmer_len_end<kmer_len) kmer_len_end=kmer_len+1;
	
	typedef Iterator<StringSet<String<Dna5> > >::Type TStringSetIterator;
	unsigned int slen=0;
	String<Dna5> Kmer;
	

	
	for(;kmer_len<=kmer_len_end;++kmer_len){//loop over all possible Kmer-length 
		for(TStringSetIterator it=begin(seq.seqs);it!=end(seq.seqs);++it){//loop over all sequences
			slen=length(value(it));//length of the current seq
			if(slen==0) continue;
			else if(slen<kmer_len) continue;
		
		
		
			for(unsigned int i=0;i<slen-kmer_len+1;++i){//loop over all Kmere in the current sequence
			
				Kmer=infix(value(it),i,i+kmer_len);
			
				if(seq.seqCounter.find(Kmer)!=seq.seqCounter.end()) continue;// if Kmer is in the Map -->nothing to do
				CountKmer(back,finderB,Kmer);
				CountKmer(seq,finder,Kmer);
				
			}
		}

	}

		
	
}

//gets the current Kmer and searches it in the index
//max(cumulated sum of the counter)=SeqsNumber --> counts number of sequences containing the motif
void CountKmer( Seq &seq, 
				Finder<Index<StringSet<String<Dna5> > > > &finder, 
				String<Dna5> &Kmer){
	
			std::vector<int> CounterV(seq.SeqsNumber+1,0);//counter for storing 1 or 0 for each Seq + the cumulated sum of the counter in the last field
			
			clear(finder);
			while(find(finder,Kmer)){//search the current Kmer in all sequences
				
			
				if(seq.c>1){//check only in foreground
					clear(seq.results);
					findIntervals(seq.intervalTrees[beginPosition(finder).i1], beginPosition(finder).i2, endPosition(finder).i2, seq.results);
				}// if results>0 --> do not count:Kmer is overlapping
				if(CounterV[beginPosition(finder).i1] == 0 && length(seq.results)==0){//count number of sequences containing the motif, not the occurrences to avoid problems with self-overlapping
					++CounterV[beginPosition(finder).i1];
					++CounterV[seq.SeqsNumber];//last Position in CounterV is cumulated sum
				}
			
			}
			seq.seqCounter[Kmer]=CounterV[seq.SeqsNumber];
			CounterV.clear();


}



void loopOverKmer(  Seq &seq,
					String<Iupac> &temp,
					String<Iupac> &Kmer,
					Iterator<String<Iupac> >::Type &tempIt,
					Finder<Index<StringSet<String<Dna5> > > > &finder,
					unsigned int &counter,
					std::vector<int> &CounterV,
					IupacMaps &IMap){
	
	String<Dna5> replace;
	Iterator<String<Dna5> >::Type replaceIt;
	Iterator<String<Iupac> >::Type tempIttemp;
	char resetTemp;

	if(tempIt==end(temp)) return;
	if((*tempIt == 'A' || *tempIt == 'C' ||*tempIt == 'G' ||*tempIt == 'T')){
		loopOverKmer(seq,temp,temp,++tempIt,finder,counter,CounterV,IMap);//only replace the position with a wildcard
		return;
	}
	replace=IMap.IupacMapReplaceReversed[*tempIt];
	
	replaceIt = begin(replace);	
			
	for(;replaceIt!=end(replace);++replaceIt){
		temp = Kmer;// reset temp
		resetTemp = IMap.IupacMapInversed[replace];
		*tempIt = *replaceIt;
		//if not end, call fkt. with temp 
		// if end, call find --> &counter
		tempIttemp=tempIt;//the  rekursive function-call with tempIttemp, because the loop needs tempIt 
		if(tempIt+1!=end(temp)){
			loopOverKmer(seq,temp,temp,++tempIttemp,finder,counter,CounterV,IMap);

		}
		
		if((*tempIttemp == 'A' || *tempIttemp == 'C' || *tempIttemp == 'G' || *tempIttemp == 'T' || tempIttemp == end(temp))){
					
			
			FindKmer(seq,temp,finder,counter,CounterV);
		}
				
		*tempIt=resetTemp;			
					
	}
				
			



}

void FindKmer(  Seq &seq,
				String<Iupac> &temp,
				Finder<Index<StringSet<String<Dna5> > > > &finder,
				unsigned int &counter,
				std::vector<int> &CounterV){



					clear(finder);
					while(find(finder,temp)){//search the current Kmer in all sequences
						
						
						if(seq.c>1){//check only in foreground
							clear(seq.results);
							findIntervals(seq.intervalTrees[beginPosition(finder).i1], beginPosition(finder).i2, endPosition(finder).i2, seq.results);
						}
						if(CounterV[beginPosition(finder).i1] == 0 && length(seq.results)==0){//count number of sequences containing the motif, not the occurrences to avoid problems with self-overlapping
							
							++CounterV[beginPosition(finder).i1];
							++CounterV[seq.SeqsNumber];//last Position in CounterV is cumulated sum
							++counter;
							
						}
			
					}


}





/*
	for counting the top 100 generalizedKmere exact
	todo -->template
*/
void CountKmer( std::map<String<Iupac>,unsigned int > &seqCounter, 
				Finder<Index<StringSet<String<Dna5> > > > &finder, 
				String<Iupac> &Kmer,
				Seq &seq,
				IupacMaps &IMap){
			
	
			
			String<Iupac> temp;	
			Iterator<String<Iupac> >::Type tempIt;
			
			temp=Kmer;
			tempIt = begin(temp);
			unsigned int counter=0;
			std::vector<int> CounterV(seq.SeqsNumber+1,0);//counter for storing 1 or 0 for each Seq + the cumulated sum of the counter in the last field
			
			
			loopOverKmer(seq,temp,Kmer,tempIt,finder,counter,CounterV,IMap);

			
			
			seqCounter[Kmer]=counter;
			CounterV.clear();
			
			/*
				Der clear hier bewirkt, dass bei ASG entweder ACG oder AGG vorkommen darf,
				falls beide vorkommen(in einer Sequenz) wird der counter trotzdem nur um 1 erhöht
				-->Counter für AGG ist falsch, ASG stimmt jedoch wieder. Counter[AGG] --> irrelevant
				-->AGG wird jedes mal wenns benötigt wird neu berechnet-->optimieren
			*/
		

}



//void replaceKmer(   Seq &seq,
//					String<unsigned int> &replaceString){
//
//	unsigned int beg =0;
//	unsigned int endP =0;
//	unsigned int Snumber =0;
//
//	//std::cout<<stringNumber<<" "<<begin<<" "<<end<<std::endl;
//	Iterator<String<unsigned int> >::Type StringIt;
//	StringIt = begin(replaceString);
//	for(;StringIt!=end(replaceString);++StringIt){
//		Snumber = *StringIt;
//		++StringIt;
//		beg		= *StringIt;
//		++StringIt;
//		endP		= *StringIt;
//		//std::cout<<Snumber<<" "<<beg<<" "<<endP<<std::endl;
//		for(;beg<endP;++beg)
//		{
//			seq.seqs[Snumber][beg]='N';
//		}
//			//replace(seq.seqs[*StringIt],*(++StringIt),*(++StringIt),'N');
//	}
//
//
//}


///////////////////////////////////
void FindTopKmer(Seq &seq,
				String<Iupac> &temp,
				Finder<Index<StringSet<String<Dna5> > > > &finder,
				unsigned int &counter,
				std::vector<int> &CounterV){

	typedef IntervalAndCargo<unsigned, unsigned> TInterval;				
	clear(finder);
	
	while(find(finder,temp)){//search the current Kmer in all sequences
		if(seq.c==1){//if c==1 save in Intervals, otherwise add direct to Tree 
			appendValue(seq.intervals[beginPosition(finder).i1], TInterval(beginPosition(finder).i2, endPosition(finder).i2, 0)); 
		}
		else if(seq.c>1)
			addInterval(seq.intervalTrees[beginPosition(finder).i1], TInterval(beginPosition(finder).i2, endPosition(finder).i2, 0));
		
		if(CounterV[beginPosition(finder).i1] == 0){//count number of sequences containing the motif, not the occurrences to avoid problems with self-overlapping
			++CounterV[beginPosition(finder).i1];
			++CounterV[seq.SeqsNumber];//last Position in CounterV is cumulated sum
			++counter;
		}
			
	}
	for( unsigned int k =0;k< length(temp);++k){//computes frequency of the nucleotides, each motif only once per sequence
		seq.freqMatrix[k][temp[k]]+=CounterV[seq.SeqsNumber]; //GCAGCA -->adds the counter of the motif to each nucleotide-counter,the same for the next motif GCAGTA --> etc. Non-Wildcards probability=1
	}
	
	CounterV.clear();
}

void FindTopKmer(Seq &seq,
				String<Dna5> &temp,
				Finder<Index<StringSet<String<Dna5> > > > &finder,
				unsigned int &counter,
				std::vector<int> &CounterV){

	typedef IntervalAndCargo<unsigned, unsigned> TInterval;				
	clear(finder);
	
	while(find(finder,temp)){//search the current Kmer in all sequences
		if(seq.c==1){//if c==1 save in Intervals, otherwise add direct to Tree 
			appendValue(seq.intervals[beginPosition(finder).i1], TInterval(beginPosition(finder).i2, endPosition(finder).i2, 0)); 
		}
		else if(seq.c>1)
			addInterval(seq.intervalTrees[beginPosition(finder).i1], TInterval(beginPosition(finder).i2, endPosition(finder).i2, 0));
		
		if(CounterV[beginPosition(finder).i1] == 0){//count number of sequences containing the motif, not the occurrences to avoid problems with self-overlapping
			++CounterV[beginPosition(finder).i1];
			++CounterV[seq.SeqsNumber];//last Position in CounterV is cumulated sum
			++counter;
		}
			
	}
	for( unsigned int k =0;k< length(temp);++k){//computes frequency of the nucleotides, each motif only once per sequence
		seq.freqMatrix[k][temp[k]]+=CounterV[seq.SeqsNumber];//GCAGCA -->adds the counter of the motif to each nucleotide-counter,the same for the next motif GCAGTA --> etc. Non-Wildcards probability=1
	}
	
	CounterV.clear();
}


void loopOverTopKmer(	Seq &seq,
						String<Iupac> &temp,
						String<Iupac> &Kmer,
						Iterator<String<Iupac> >::Type &tempIt,
						Finder<Index<StringSet<String<Dna5> > > > &finder,
						unsigned int &counter,
						std::vector<int> &CounterV,
						IupacMaps &IMaps){

						
			String<Dna5> replace;
			Iterator<String<Dna5> >::Type replaceIt;
			Iterator<String<Iupac> >::Type tempIttemp;
			
			char resetTemp;//bei mehr als einer wildcard, müssen die weiter hinten liegenden nach Abarbeitung resetet werden, ansonsten werden diese im nächsten schritt übergangen

			if(tempIt==end(temp)) return;
			
			if((*tempIt == 'A' || *tempIt == 'C' ||*tempIt == 'G' ||*tempIt == 'T')){
				loopOverTopKmer(seq,temp,temp,++tempIt,finder,counter,CounterV,IMaps);//only replace the position with a wildcard
				return;//always return-->loop otherwise
			}
				
			
				replace=IMaps.IupacMapReplaceReversed[*tempIt];
				replaceIt = begin(replace);	
				for(;replaceIt!=end(replace);++replaceIt){
					
					temp = Kmer;// reset temp
					resetTemp = IMaps.IupacMapInversed[replace]; //if Y gets replaced-->  replace = CT -->  resetTemp = Y
					
					*tempIt = *replaceIt;
					
					tempIttemp=tempIt;//the rekursive function call with tempIttemp--> loop needs tempIt 
			
					
					if(tempIt+1!=end(temp)){
						
						loopOverTopKmer(seq,temp,temp,++tempIttemp,finder,counter,CounterV,IMaps);
						
					}
					
					if((*tempIttemp == 'A' || *tempIttemp == 'C' || *tempIttemp == 'G' || *tempIttemp == 'T' || tempIttemp == end(temp))){
						/**
						posTemp gibt die Stelle an, an der man sich grad befindet 
						--> freqMatrix enthält für jede Position des Motivs die W'keit für ACGT
						--> counter ist die Anzahl wie oft das Wildcard Motiv insgesamt gefunden wurde
						--> CounterV an letzter Stelle die Anzahl für das jeweilige nicht-Wildcard Motiv
						--> bei GSAKYA z.B. als Motiv wird jedes Motiv bei 'S' vier mal gesucht(durch die anderen 2 Wildcards K und Y)
						--> CounterV für 1 bis posTemp aufaddieren --> in freqMatrix und zwar für die jeweiligen *tempIt-chars
						--> am Ende alle durch counter teilen --> aufpassen, für jeweilige pos gibts verschiedene counter
						--> FindKmer wird nur mit ganzen aufgerufen, also alle addieren, dann ist der counter auch gleich?
						**/
						std::vector<int> CounterV(seq.SeqsNumber+1,0);
						FindTopKmer(seq,temp,finder,counter,CounterV);
						
					}
					*tempIt=resetTemp;
					
				}


}

/***
	Computes PWM 
***/
void BuildFrequencyMatrix(  Finder<Index<StringSet<String<Dna5> > > > &finder,
							String<Iupac> &Kmer,
							Seq &seq, 
							IupacMaps &IMaps){
			
			//freqMatrix -->unsigned int = position in Kmer; position 1 in map = prob. for A, pos. 2 = prob. for C...
			appendValue(seq.allMotifs,Kmer,Exact());
								
			String<Iupac> temp;	
			
			Iterator<String<Iupac> >::Type tempIt;
			
			temp=Kmer;
			tempIt = begin(temp);
			unsigned int counter=0;
			
			std::vector<int> CounterV(seq.SeqsNumber+1,0);
		
			
			loopOverTopKmer(seq,temp,Kmer,tempIt,finder,counter,CounterV,IMaps);
			CounterV.clear();
			
			if(counter>0){//normalize counter
				
				
				for( unsigned int k =0;k< length(temp);++k){
					
					seq.freqMatrix[k]['A']=(seq.freqMatrix[k]['A']+seq.frequencies[0])/(counter+1);//corrected freq (with pseudocount)
					seq.freqMatrix[k]['C']=(seq.freqMatrix[k]['C']+seq.frequencies[1])/(counter+1); 
					seq.freqMatrix[k]['G']=(seq.freqMatrix[k]['G']+seq.frequencies[2])/(counter+1); 
					seq.freqMatrix[k]['T']=(seq.freqMatrix[k]['T']+seq.frequencies[3])/(counter+1);  
				}
	

				}
			else
				seq.freqMatrix.clear();
				
				
				
			
}

/**
		Computes PWM --> not-generalized Kmer-Version
**/
void BuildFrequencyMatrix(  Finder<Index<StringSet<String<Dna5> > > > &finder,
							String<Dna5> &Kmer,
							Seq &seq, 
							IupacMaps &IMaps){
			
			//freqMatrix -->unsigned int = position in Kmer, position 1 in map = prob. for A, pos. 2 = prob. for C...
			
			appendValue(seq.allMotifs,Kmer,Exact());
								
			String<Iupac> temp;	
			Iterator<String<Iupac> >::Type tempIt;


			
			temp=Kmer;
			tempIt = begin(temp);
			unsigned int counter=0;
			
			std::vector<int> CounterV(seq.SeqsNumber+1,0);
			FindTopKmer(seq,Kmer,finder,counter,CounterV);
			CounterV.clear();
			
			
				
			if(counter>0){//normalize counter
				
				
				for( unsigned int k =0;k< length(temp);++k){
					
					seq.freqMatrix[k]['A']=(seq.freqMatrix[k]['A']+seq.frequencies[0])/(counter+1);//corrected freq (with pseudocount)
					seq.freqMatrix[k]['C']=(seq.freqMatrix[k]['C']+seq.frequencies[1])/(counter+1); 
					seq.freqMatrix[k]['G']=(seq.freqMatrix[k]['G']+seq.frequencies[2])/(counter+1); 
					seq.freqMatrix[k]['T']=(seq.freqMatrix[k]['T']+seq.frequencies[3])/(counter+1); 
				}
	

				}
			else
				seq.freqMatrix.clear();
				
				
				
			
}
//////////////////////////////////////

//FreqMatrix output
void PrintMap(  Seq &seq,
				bool foreground){

	std::map<Iupac,double> freq;
	std::cout<<std::endl;
	if(foreground)
		std::cout<<"foreground: "<<std::endl;
	else
		std::cout<<"background: "<<std::endl;
	for(unsigned int j=0;j<length(seq.freqMatrix);++j){
				freq=seq.freqMatrix[j];
				std::cout<<"Position: "<<j<<" A: "<<freq['A']<<std::endl;
				std::cout<<"Position: "<<j<<" C: "<<freq['C']<<std::endl;
				std::cout<<"Position: "<<j<<" G: "<<freq['G']<<std::endl;
				std::cout<<"Position: "<<j<<" T: "<<freq['T']<<std::endl;
				std::cout<<std::endl;
			}
			

}



void saveData(	Seq &seq,
				std::ofstream &PWM,
				unsigned c){

	

	String<char> pwm;
	char *buffer= new char[33];
	append(pwm,"C:/Users/David/Desktop/PWM");
	//append(pwm,sprintf(buffer,"%d",c));
	append(pwm,itoa(c,buffer,10));
	append(pwm,".txt");
	
	
	PWM.open(toCString(pwm),std::ios::out|std::ios::app);
	
	for(unsigned i=0;i<length(seq.freqMatrix);++i){
		PWM<<seq.freqMatrix[i]['A'];
		write(PWM," ");
		

	}
	write(PWM,"\r\n");// windows
	PWM.close();
	PWM.open(toCString(pwm),std::ios::out|std::ios::app);
	
	for(unsigned i=0;i<length(seq.freqMatrix);++i){
		PWM<<seq.freqMatrix[i]['C'];
		write(PWM," ");

	}
	write(PWM,"\r\n");
	PWM.close();
	PWM.open(toCString(pwm),std::ios::out|std::ios::app);
	
	for(unsigned i=0;i<length(seq.freqMatrix);++i){
		PWM<<seq.freqMatrix[i]['G'];
		write(PWM," ");

	}
	write(PWM,"\r\n");
	PWM.close();
	PWM.open(toCString(pwm),std::ios::out|std::ios::app);
	
	for(unsigned i=0;i<length(seq.freqMatrix);++i){
		PWM<<seq.freqMatrix[i]['T'];
		write(PWM," ");

	}
	write(PWM,"\r\n");
	write(PWM,"\r\n");
	PWM.close();
	delete[] buffer;

}

/*****
		Computes weightMatrix -->needed for InformationContentMatrix
*****/
void BuildWeightMatrix(Seq &seq){

	
	
	for(unsigned int j=0;j<length(seq.freqMatrix);++j){
		seq.weightMatrix[j]['A'] = log(seq.freqMatrix[j]['A']/seq.frequencies[0]);
		seq.weightMatrix[j]['C'] = log(seq.freqMatrix[j]['C']/seq.frequencies[1]);
		seq.weightMatrix[j]['G'] = log(seq.freqMatrix[j]['G']/seq.frequencies[2]);
		seq.weightMatrix[j]['T'] = log(seq.freqMatrix[j]['T']/seq.frequencies[3]);
	}




}
/*****
		computes the information content matrix of a motif
*****/

void BuildInformationContentMatrix(Seq &seq){

	double e = 3/(2*log(2)*seq.SeqsNumber);//4-1 = 3, 4= Anzahl Nukleotide --> small-error correction

	for(unsigned int j=0;j<length(seq.freqMatrix);++j){
		seq.InformationContentMatrix[j]['A'] = seq.freqMatrix[j]['A']*seq.weightMatrix[j]['A'];
		seq.InformationContentMatrix[j]['C'] = seq.freqMatrix[j]['C']*seq.weightMatrix[j]['C'];
		seq.InformationContentMatrix[j]['G'] = seq.freqMatrix[j]['G']*seq.weightMatrix[j]['G'];
		seq.InformationContentMatrix[j]['T'] = seq.freqMatrix[j]['T']*seq.weightMatrix[j]['T'];
		seq.seqLogoMatrix[j]= 2- (-(seq.InformationContentMatrix[j]['A']+seq.InformationContentMatrix[j]['C']+seq.InformationContentMatrix[j]['G']+seq.InformationContentMatrix[j]['T'])+ e);
	}

}

double ComparePWM(	std::map<Iupac,double>  &freqMatrix1,
					std::map<Iupac,double>  &freqMatrix2, 
					Entropy const & tag){


	
	/****
			Für die übergebenen Spalten der Matrizen wird die Entropy(=Kullback-Leibler-Divergenz) berechnet
			Eintrag ist = 0, wenn die Werte identisch sind. Je größer die Zahl, desto unterschiedlicher die Werte
			Was wenn Wert = Hintergrundverteilung? --> Eintrag wäre = 0, obwohl das nichts mit dem Motiv zu tun hat

			Für jede Spalte berechnen und im Nachhinein durch Spaltenanzahl teilen 
	****/
	double columnEntropy=0;
	
	
	columnEntropy  = freqMatrix1['A']*log(freqMatrix1['A']/freqMatrix2['A']);
	columnEntropy += freqMatrix1['C']*log(freqMatrix1['C']/freqMatrix2['C']);
	columnEntropy += freqMatrix1['G']*log(freqMatrix1['G']/freqMatrix2['G']);
	columnEntropy += freqMatrix1['T']*log(freqMatrix1['T']/freqMatrix2['T']);
	columnEntropy += freqMatrix2['A']*log(freqMatrix2['A']/freqMatrix1['A']);
	columnEntropy += freqMatrix2['C']*log(freqMatrix2['C']/freqMatrix1['C']);
	columnEntropy += freqMatrix2['G']*log(freqMatrix2['G']/freqMatrix1['G']);
	columnEntropy += freqMatrix2['T']*log(freqMatrix2['T']/freqMatrix1['T']);
	columnEntropy  = columnEntropy/2;

	
	return columnEntropy;

}
/****
			Computes the euclidean distance
		
****/
double ComparePWM(	std::map<Iupac,double>  &freqMatrix1,
					std::map<Iupac,double>  &freqMatrix2, 
					Euclid const & tag){
	

	double columnEntropy = 0;

	columnEntropy  = (freqMatrix1['A'] - freqMatrix2['A'])*(freqMatrix1['A'] - freqMatrix2['A']);
	columnEntropy += (freqMatrix1['C'] - freqMatrix2['C'])*(freqMatrix1['C'] - freqMatrix2['C']);
	columnEntropy += (freqMatrix1['G'] - freqMatrix2['G'])*(freqMatrix1['G'] - freqMatrix2['G']);
	columnEntropy += (freqMatrix1['T'] - freqMatrix2['T'])*(freqMatrix1['T'] - freqMatrix2['T']);
	columnEntropy  = sqrt(columnEntropy);
	columnEntropy  = columnEntropy/4;
	
	return columnEntropy;

}
/****
		computes the Pearson correlation coefficient
****/
double ComparePWM(	std::map<Iupac,double>  &freqMatrix1,
					std::map<Iupac,double>  &freqMatrix2, 
					Pearson const & tag){

	double columnEntropy = 0;
	double columnEntropyDivisorA = 0;
	double columnEntropyDivisorB = 0;
	double Mean1 = 0;
	double Mean2 = 0;

	Mean1= (freqMatrix1['A'] + freqMatrix1['C'] + freqMatrix1['G'] + freqMatrix1['T'])/4;
	Mean2= (freqMatrix2['A'] + freqMatrix2['C'] + freqMatrix2['G'] + freqMatrix2['T'])/4;

	columnEntropy += (freqMatrix1['A'] - Mean1)*(freqMatrix2['A'] - Mean2);
	columnEntropy += (freqMatrix1['C'] - Mean1)*(freqMatrix2['C'] - Mean2);
	columnEntropy += (freqMatrix1['G'] - Mean1)*(freqMatrix2['G'] - Mean2);
	columnEntropy += (freqMatrix1['T'] - Mean1)*(freqMatrix2['T'] - Mean2);

	columnEntropyDivisorA += (freqMatrix1['A'] - Mean1)*(freqMatrix1['A'] - Mean1);
	columnEntropyDivisorA += (freqMatrix1['C'] - Mean1)*(freqMatrix1['C'] - Mean1);
	columnEntropyDivisorA += (freqMatrix1['G'] - Mean1)*(freqMatrix1['G'] - Mean1);
	columnEntropyDivisorA += (freqMatrix1['T'] - Mean1)*(freqMatrix1['T'] - Mean1);

	columnEntropyDivisorB += (freqMatrix2['A'] - Mean2)*(freqMatrix2['A'] - Mean2);
	columnEntropyDivisorB += (freqMatrix2['C'] - Mean2)*(freqMatrix2['C'] - Mean2);
	columnEntropyDivisorB += (freqMatrix2['G'] - Mean2)*(freqMatrix2['G'] - Mean2);
	columnEntropyDivisorB += (freqMatrix2['T'] - Mean2)*(freqMatrix2['T'] - Mean2);

	columnEntropyDivisorA = columnEntropyDivisorA * columnEntropyDivisorB;
	columnEntropyDivisorA = sqrt(columnEntropyDivisorA);

	columnEntropy = columnEntropy/columnEntropyDivisorA;


	return columnEntropy;

}

/****
		Computes the mean of the two matrices 
****/
void BuildMeanOf2PWMs(	Seq &seq,
						std::map<unsigned int,std::map<Iupac,double> > &freqMatrix1,
						std::map<unsigned int,std::map<Iupac,double> > &freqMatrix2){

	for(unsigned int j=0;j<length(freqMatrix1);++j){
		freqMatrix1[j]['A'] = (freqMatrix1[j]['A']+freqMatrix2[j]['A'])/2;
		freqMatrix1[j]['C'] = (freqMatrix1[j]['C']+freqMatrix2[j]['C'])/2;
		freqMatrix1[j]['G'] = (freqMatrix1[j]['G']+freqMatrix2[j]['G'])/2;
		freqMatrix1[j]['T'] = (freqMatrix1[j]['T']+freqMatrix2[j]['T'])/2;
	}
	

}
/****
			Updates the distance with Complete Linkage
			Replaces PWM x with the new one
			
****/

void UpdateDistantMatrix(	int n,
							int x, 
							int y, 
							std::vector<std::vector<String<double> > > &compare, 
							CompleteLinkage const & tag ){

	
	unsigned j;
	for (j = 0; j < x; j++)
		compare[x][j][0] = std::max(compare[y][j][0],compare[x][j][0]);
	for (j = x+1; j < y; j++)
		compare[j][x][0] = std::max(compare[y][j][0],compare[j][x][0]);
	for (j = y+1; j < n; j++)
		compare[j][x][0] = std::max(compare[j][y][0],compare[j][x][0]);
	for (j = 0; j < y; j++) compare[y][j][0] = compare[n-1][j][0];
	for (j = y+1; j < n-1; j++) compare[j][y][0] = compare[n-1][j][0];




}
/****
			Updates the distance with Single Linkage
			Replaces PWM x with the new one
			
****/
void UpdateDistantMatrix(	int n,
							int x, 
							int y, 
							std::vector<std::vector<String<double> > > &compare, 
							SingleLinkage const & tag ){

	
	unsigned j;
	for (j = 0; j < x; j++)
		compare[x][j][0] = std::min(compare[y][j][0],compare[x][j][0]);
	for (j = x+1; j < y; j++)
		compare[j][x][0] = std::min(compare[y][j][0],compare[j][x][0]);
	for (j = y+1; j < n; j++)
		compare[j][x][0] = std::min(compare[j][y][0],compare[j][x][0]);
	for (j = 0; j < y; j++) compare[y][j][0] = compare[n-1][j][0];
	for (j = y+1; j < n-1; j++) compare[j][y][0] = compare[n-1][j][0];




}


/****
		Updates the distance with Average Linkage
****/

void UpdateDistantMatrix(	int n, 
							int x, 
							int y,
							std::vector<std::vector<String<double> > > &compare,
							std::vector<unsigned> &weights, 
							AverageLinkage const & tag ){


	/***
			Die alten Einträge gewichtet addieren --> Gewichte am Anfang =1 ^= Clustergröße
			Zur Normierung durch die addierten Gewichte teilen
	***/
	unsigned j;
	unsigned sumOfweights = weights[x] + weights[y];
	
	for (j = 0; j < x; j++)
		compare[x][j][0] = (compare[y][j][0]*weights[y]+compare[x][j][0]*weights[x])/sumOfweights;
	for (j = x+1; j < y; j++)
		compare[j][x][0] = (compare[y][j][0]*weights[y]+compare[j][x][0]*weights[x])/sumOfweights;
	for (j = y+1; j < n; j++)
		compare[j][x][0] = (compare[j][y][0]*weights[y]+compare[j][x][0]*weights[x])/sumOfweights;
	for (j = 0; j < y; j++) compare[y][j][0] = compare[n-1][j][0];
	for (j = y+1; j < n-1; j++) compare[j][y][0] = compare[n-1][j][0];
	weights[x]=sumOfweights;
	for(j=0;j+y<length(weights)-1;++j)
			weights[y+j]=weights[y+j+1];

}

/****
		computes the minimum within the matrix 'compare'
****/
void minDifferenceInMatrix(	unsigned n,
							String<double> &minDifference,
							std::vector<std::vector<String<double> > > compare){

	minDifference[0]=0;
	for(unsigned i=0;i<n-1;++i){

		for(unsigned j=i+1;j<n;++j){
			
			if(minDifference[0]==0 || minDifference[0]>compare[j][i][0]){
				minDifference[0]=compare[j][i][0];
				minDifference[1]=i;//PWM i
				minDifference[2]=j;//PWM j --> PWM i and j = PWMs of the smallest value in the matrix
				
			}

		}

	}

}

/****
		computes the distant matrix with the alignment of the PWMs --> for overlapping PWMs
****/
void computesDistantMatrix(	String< std::map<unsigned int,std::map<Iupac,double> > > &allPWMs,
							std::vector<std::vector<String<double> > > &compare, 
							unsigned allPWMsLength){

	unsigned j;
	unsigned i;
	for(i=0;i<allPWMsLength-1;++i){//computes distances with ComparePWM and saves it in the matrix compare

		for( j=i+1;j<allPWMsLength;++j){

			
			compare[j][i]=AlignPWMs(allPWMs[i],allPWMs[j]);
			
			

		}
		

	}

}
/****
		Dr= the sum of all distances for all points in cluster r
****/

double computeDr(	Cluster &cluster,
					std::vector<std::vector<String<double> > > copyCompare){

	double Dr=0;
	
	for(unsigned i=0;i<length(cluster.content);++i){
		for(unsigned j=0;j<length(cluster.content);++j){
			if(cluster.content[i]<cluster.content[j]){
				
				Dr+= copyCompare[cluster.content[j]][cluster.content[i]][0];//copyCompare[i][j] --> i>j, therefore reverse if i<j 
				
			}
			else if(cluster.content[j]<cluster.content[i]){
				
				Dr+= copyCompare[cluster.content[i]][cluster.content[j]][0];
				
			}
		}
	}
	
	return Dr/2;
}
/****
		Wk = the pooled within-cluster sum of squares around the cluster means
****/

double computeWk(	int n,
					String<Cluster> &cluster){// n=allPWMsLength-n

	/***
		Wenn ein neues Cluster hinzukommt einfach addieren/n
		Wenn ein bestehendes verändert wird, den Wert austauschen
	***/


	if(n>0){
		cluster[n].Wk=cluster[n-1].Wk;
	}
	else
		cluster[n].Wk=0;

	if(cluster[n].left>=0 && cluster[n].right>=0) 
		cluster[n].Wk +=cluster[n].Dr/length(cluster[n].content);//no need to divide by 2 --> content contains only the half
	else{
		if(cluster[n].left<0)
			cluster[n].Wk=cluster[n].Wk - cluster[-cluster[n].left-1].Dr/length(cluster[-cluster[n].left-1].content);
		if(cluster[n].right<0)
			cluster[n].Wk=cluster[n].Wk - cluster[-cluster[n].right-1].Dr/length(cluster[-cluster[n].right-1].content);	

		cluster[n].Wk= cluster[n].Wk + cluster[n].Dr/length(cluster[n].content);
	}
	
	return cluster[n].Wk;


}


/***
	Erzeugt gleichverteilte Daten, die zur Berechnung der Gap herangezogen werden
	Wird B mal aufgerufen für verschiedene Data-Sets

***/
void computeReferenceData(String< std::map<unsigned int,std::map<Iupac,double> > > &Reference,
						  unsigned allPWMsLength,
						  unsigned PWMLength ){

	/***
			Berechnet für jede Spalte eines PWMs gleichverteilte Zufallszahlen
			Erzeugt so viele PWMs wie auch gemessen wurden

			Spalten:
			Zufällig zwischen ACGT auswählen --> x von 0 bis 0.99
			Nächste aus ACGT --> y von 0 bis (1-x)
			Nächstes --> z von 0 bis (1-x-y)

			PWMLength = Länge der Motive

	***/
	
	reserve(Reference,allPWMsLength);
	std::map<unsigned int,std::map<Iupac,double> > ReferenzFreq;

	String<Dna> ACGT;
	int pos;
	const int SEED = 7;
	Rng<MersenneTwister> rng(SEED);
	

	double uniDouble;
	

	for(unsigned i=0;i<allPWMsLength;++i){//genauso viele Reference-Daten erstellen, wie PWMs vorhanden
		for(unsigned j=0;j<PWMLength;++j){
			ACGT="ACGT";
			uniDouble=0;
			for(unsigned k=0;k<3;++k){
				
				Pdf<Uniform<int> > uniformInt(0, length(ACGT)-1);
			
				Pdf<Uniform<double> > uniformDouble(0, 0.99-uniDouble);
				pos=pickRandomNumber(rng, uniformInt);
				ReferenzFreq[j][ACGT[pos]]=pickRandomNumber(rng, uniformDouble);
				uniDouble+=ReferenzFreq[j][ACGT[pos]];
				erase(ACGT,pos);
			}
			ReferenzFreq[j][ACGT[0]]=1-uniDouble;
			clear(ACGT);

		} 
		
		appendValue(Reference,ReferenzFreq,Exact());
		
		clear(ReferenzFreq);
	}



}
/****
		l_quer= the mean of the sum of all Wks from the random reference data
****/
void compute_l_quer(String<Cluster> &cluster, 
					String<double> &l_quer_for_k,
					String<String<double> > &allWk,
					unsigned i,
					int B,
					unsigned allPWMsLength){

		unsigned k=0;
		
		for(unsigned j=0;j<allPWMsLength-1;++j){//adds the log(Wk's) from all B-Reference-Data -->only the sum is required
			/*****
					j != number of cluster -->allPWMsLength - (j+1) =number of Cluster
					
			*****/
			k=allPWMsLength-(j+1);
			appendValue(allWk[i],cluster[j].Wk);
			
			l_quer_for_k[0]+=log(cluster[j].Wk)/B;

			if(i==0){
				l_quer_for_k[k]=log(cluster[j].Wk)/B;
			}
			else{
				l_quer_for_k[k]+=log(cluster[j].Wk)/B;
			}
		}

		

}

void computeGap(String<double> &l_quer_for_k,
				String<Cluster> &observedCluster,
				String<double> &Gap,
				unsigned allPWMsLength){

	
	unsigned k=0;
	for(unsigned j=0;j<allPWMsLength-1;++j){
		k=allPWMsLength-(j+1);
		Gap[k]=l_quer_for_k[k]-log(observedCluster[j].Wk);
	}

}

/****
		computes the standard deviation sdk and sk= sdk*sqrt(1 + 1/B)
****/
void compute_sdk_and_sk(String<double> &sk,
						String<double> &l_quer_for_k,
						String<String<double> > &allWk,
						unsigned allPWMsLength,
						int B){

	String<double> sdk;//standard deviation
	resize(sdk,allPWMsLength);
	unsigned k=0;
	for(unsigned i=0;i<B;++i){

		for(unsigned j=0;j<allPWMsLength-1;++j){
			k=allPWMsLength-(j+1);
			if(i==0)
				sdk[k]=(log(allWk[i][j])-l_quer_for_k[0])*(log(allWk[i][j])-l_quer_for_k[0]);
			else
				sdk[k]+=(log(allWk[i][j])-l_quer_for_k[0])*(log(allWk[i][j])-l_quer_for_k[0]);
		}
	}

	
	

	/*******************
		compute sk
	********************/

	double x=1+1/double(B);
	for(unsigned j=0;j<allPWMsLength-1;++j){
		k=allPWMsLength-(j+1);
		sdk[k]=sqrt(sdk[k]/B);
		sk[k]=sdk[k]*sqrt(x);

	}
	clear(sdk);

}

unsigned computeGapStat(String< std::map<unsigned int,std::map<Iupac,double> > > &allPWMs,
						String<Cluster> &observedCluster){

	String<Cluster> cluster;
	String<String<double> > allWk;
	String<double> Gap;
	String<double> l_quer_for_k;
	String<int> traceback;
	String<double> sk;
	int B=5;//B reference-data sets 


	String< std::map<unsigned int,std::map<Iupac,double> > > Reference;
	unsigned allPWMsLength=length(allPWMs);
	unsigned PWMLength=length(allPWMs[0]);
	
	resize(allWk,B);
	
	resize(sk,allPWMsLength);
	resize(Gap,allPWMsLength);
	resize(l_quer_for_k,allPWMsLength);
	l_quer_for_k[0]=0;
	
	/************
		compute l_quer 
	************/
	for(unsigned i=0;i<B;++i){

		computeReferenceData( Reference,allPWMsLength, PWMLength );
		
		PWMClustering(Reference,cluster,traceback);
		compute_l_quer(cluster,l_quer_for_k,allWk,i,B,allPWMsLength);
		clear(Reference);
		clear(cluster);
	}
	
	/******************
		compute Gap 
	******************/
	computeGap(l_quer_for_k,observedCluster,Gap,allPWMsLength);


	/******************
		compute sdk && sk 
	******************/

	compute_sdk_and_sk(sk,l_quer_for_k,allWk,allPWMsLength,B);

	/*******************
		compute smallest k 
		--> returns optimal cluster-size
	********************/
	unsigned j=0;
	for(unsigned k=1;k<allPWMsLength-1;++k){
		j=allPWMsLength-(k+1);
		if(Gap[j]<=Gap[j+1]-sk[j+1]){
			std::cout<<"j "<<j<<std::endl;
			std::cout<<j+1<<" Cluster sind das Optimum, also nach dem "<<allPWMsLength - (j + 1)<<" Cluster-Schritt"<<std::endl;
			clear(allWk);
			clear(Gap);
			return (allPWMsLength - (j + 1));
		}

	}


	
	clear(allWk);
	clear(Gap);
	clear(l_quer_for_k);
	clear(sk);
	clear(traceback);

	return 0;

}

/****
		The traceback of the clustering
		--> prints the clustered motifs in the right order
****/
void clusterTraceback(	String<String<Iupac> > &allMotifs,
						String<int> &traceback, 
						unsigned j){

	Iterator<String<int> >::Type tracebackIt;
	if(j==0){	
		std::cout<<"Kein Clustern sinnvoll"<<std::endl;
		std::cout<<allMotifs[*begin(traceback)]<<" "<<allMotifs[*++begin(traceback)]<<std::endl;
		return;
	}
	std::cout<<"Traceback: "<<std::endl;
	unsigned o=0;
	for(tracebackIt=begin(traceback);tracebackIt!=end(traceback) && o<j;++tracebackIt, ++o){

		std::cout<<allMotifs[*tracebackIt]<<" ";
		++tracebackIt;
		std::cout<<allMotifs[*tracebackIt]<<std::endl;
	}
	
}
/****
		performes the actual clustering with average linkage (currently)
****/
void PWMClustering(String< std::map<unsigned int,std::map<Iupac,double> > > &allPWMs,
				   String<Cluster> &cluster,
				   String<int> &traceback){
	
	unsigned allPWMsLength=length(allPWMs);	
	String<double> minDifference;
	
	std::vector<std::vector<String<double> > > compare(allPWMsLength, std::vector<String<double> >(allPWMsLength));
	std::vector<int> clusterId(allPWMsLength);

	resize(cluster,allPWMsLength);
	resize(minDifference,3);
	minDifference[0]=0;
	 
	
	unsigned j;
	for (j = 0; j < allPWMsLength; j++) clusterId[j] = j;//to assign which PWM is in which cluster
	/***
			 required by average linkage:
	***/
	std::vector<unsigned> weights(allPWMsLength,1);
	
	


	/****
			generates the distant matrix 'compare'
	****/
	computesDistantMatrix(allPWMs,compare,allPWMsLength);

	std::vector<std::vector<String<double> > > copyCompare(compare);
	String<int> temp;
	minDifferenceInMatrix(allPWMsLength,minDifference,compare);
	/****
			computes all possible cluster-sizes
	****/
	for(unsigned n=allPWMsLength;n>1;--n){
		
		
		//BuildMeanOf2PWMs(seq,seq.allPWMs[int(minDifference[1])],seq.allPWMs[int(minDifference[2])]);//Computes the mean of the two matrices
		
		UpdateDistantMatrix(n,int(minDifference[1]),int(minDifference[2]),compare,weights,AverageLinkage());//hardcoded -->ToDo
	
		
		/****
				generates vector with all objects in the cluster
				-> needed in computeDr
		****/
		if(clusterId[int(minDifference[1])]<0) 
			append(temp,cluster[-clusterId[int(minDifference[1])] - 1].content,Exact()); 
			
		else
			appendValue(temp,clusterId[int(minDifference[1])],Exact());

		if(clusterId[int(minDifference[2])]<0)
			append(temp,cluster[-clusterId[int(minDifference[2])] - 1].content,Exact()); 
		else
			appendValue(temp,clusterId[int(minDifference[2])],Exact());
			
		
		cluster[allPWMsLength-n].content=temp;
		clear(temp);
		/****
					computes Dr
		****/
		cluster[allPWMsLength-n].left=clusterId[int(minDifference[1])];
		cluster[allPWMsLength-n].right=clusterId[int(minDifference[2])];
		cluster[allPWMsLength-n].Dr=computeDr(cluster[allPWMsLength-n],copyCompare);
		
		/****
					computes Wk
		****/
		
		computeWk(allPWMsLength-n,cluster);
		

		appendValue(traceback,clusterId[int(minDifference[1])]);
		appendValue(traceback,clusterId[int(minDifference[2])]);
		/****
			change Id's 
		****/
		clusterId[int(minDifference[1])]=n-allPWMsLength-1;
		for(j=0;j+minDifference[2]<allPWMsLength-1;++j)
			clusterId[int(minDifference[2])+j]=clusterId[int(minDifference[2])+j+1];
		/****
				compute new minDifference
		****/
		minDifferenceInMatrix(n-1,minDifference,compare);
	}
	
	clear(clusterId);
	clear(minDifference);
	clear(compare);
	clear(copyCompare);
	clear(weights);

	

}


/*****Prints the Mapping:
Kmer	Seq1	Seq2	...	Seqn	CumulatedCounter
-->Template

*****/
void PrintMap(  std::map<String<Dna5>,unsigned int> &Dna5CounterMap,
				unsigned int SeqsNumber){
	std::cout<<std::endl;
	std::map<String<Dna5>,unsigned int>::iterator MapIterator;
	for(MapIterator=Dna5CounterMap.begin(); MapIterator !=Dna5CounterMap.end();++MapIterator){
		std::cout<<(*MapIterator).first<<"   ";
		
			std::cout<<(*MapIterator).second<<" ";
		std::cout<<std::endl;
	}
}

void PrintMap(  std::map<String<Iupac>,unsigned int > &Dna5CounterMap,
				unsigned int SeqsNumber){
	std::cout<<std::endl;
	std::map<String<Iupac>,unsigned int >::iterator MapIterator;
	for(MapIterator=Dna5CounterMap.begin(); MapIterator !=Dna5CounterMap.end();++MapIterator){
		std::cout<<(*MapIterator).first<<"   ";
		
			std::cout<<(*MapIterator).second<<" ";
		std::cout<<std::endl;
	}
}

void PrintMap(std::multimap<double,String<Dna5> > &pValueMap){
	std::multimap<double,String<Dna5> >::iterator MapIterator;
	for(MapIterator=pValueMap.begin();MapIterator !=pValueMap.end();++MapIterator){
		std::cout<<(*MapIterator).first<<" ";
		std::cout<<(*MapIterator).second<<std::endl;
	}
}

void PrintMap(std::multimap<double,String<Iupac> > &pValueMap){
	std::multimap<double,String<Iupac> >::iterator MapIterator;
	
	int i=0;
	for(MapIterator=pValueMap.begin();MapIterator !=pValueMap.end() && i<10;++MapIterator,++i){
		std::cout<<(*MapIterator).first<<" ";
		std::cout<<(*MapIterator).second<<std::endl;
	}
}

//to print the generalizedKmerMap
void PrintMap(std::map<String<Iupac>,unsigned int> &generalizedKmer){
	std::map<String<Iupac>, unsigned int>::iterator MapIterator;
	int i=0;
	for(MapIterator=generalizedKmer.begin();MapIterator!=generalizedKmer.end() && i<20;++MapIterator,++i){
		std::cout<<(*MapIterator).first<<" ";
		std::cout<<(*MapIterator).second<<"    ";
	}
}


////Test the Map-lengths match eachother and with the sequences
//void DebugMap(  Seq &seq,
//				Seq &back,
//				std::map<String<Dna5>,std::vector<int> > &sequencesCounter,
//				std::map<String<Dna5>,std::vector<int> > &backgroundCounter){
//
//	typedef std::map<String<Dna5>,std::vector<int> > Dna5CounterMap;
//	Dna5CounterMap::iterator MapIterator;
//	MapIterator=sequencesCounter.begin();
//	Dna5CounterMap::iterator MapIteratorB;
//	MapIteratorB=backgroundCounter.begin();
//
//	
//	SEQAN_ASSERT_EQ(length(sequencesCounter),length(backgroundCounter));
//	SEQAN_ASSERT_EQ(length((*MapIterator).second),(length(seq.ids)+1));//+1, because of the last field in vector 
//	SEQAN_ASSERT_EQ(length((*MapIteratorB).second),(length(back.ids)+1));
//
//	//std::cout<<length(sequencesCounter)<<std::endl;
//	//std::cout<<length((*MapIterator).second)<<std::endl;
//	//std::cout<<length(backgroundCounter)<<std::endl;
//	//std::cout<<length((*MapIteratorB).second)<<std::endl;
//	//std::cout<<length(seq.ids)<<std::endl;
//	//std::cout<<length(back.ids)<<std::endl;
//}
//
//void DebugMultiMap( std::map<String<Dna5>,std::vector<int> > &sequencesCounter,
//					std::multimap<double,String<Dna5> > &SortedPValue){
//	SEQAN_ASSERT_EQ(length(sequencesCounter),SortedPValue.size());
//
//}

/****
		computes all needed logFactorial-values and saves them in p
		only computed once
****/
void logFactorial(unsigned int len){
	double* p;
	unsigned int i=1;
	p = (double*)malloc(sizeof(double)*(len+1));
	p[0]=0;
		for(;i<=len;++i){
			p[i]= p[i-1] + log(i);
		}

	

	logFac =p;
}
/****
		computes the probability with the  hypergeometric distribution --> required by the fisher-exact-test 
****/
double calcFET( unsigned int a,
				unsigned int b,
				unsigned int c,
				unsigned int d){

	return exp(logFac[a+b] + logFac[c+d] + logFac[a+c] + logFac[b+d] - 
		(logFac[a+b+c+d] + logFac[a]+ logFac[b] + logFac[c] +logFac[d]));
}

/****
		modifys  a,b,c and d to be more extreme--> needed for the one-sided test
****/
void modifyFET( unsigned int a,
				unsigned int b,
				unsigned int c,
				unsigned int d, 
				double &pValue){

	
	
	pValue= calcFET(a,b,c,d);
	
		while(b!=0 && c!=0){//modify to be more extrem
			++a;
			--b;
			--c;
			++d;
			pValue += calcFET(a,b,c,d);
	
		}

}




/***************************

log((a+b)!(c+d)!(a+c)!(b+d)!/a!b!c!d!n!) = logFactorial(a+b) + logFactorial(c+d) +
										   logFactorial(a+c) + logFactorial(b+d) -
										   (logFactorial(a+b+c+d) + logFactorial(a)+
										   logFactorial(b) + logFactorial(c) +
										   logFactorial(d))

pValue = exp(log((a+b)!(c+d)!(a+c)!(b+d)!/a!b!c!d!n!))

a= Found in Sequence 	b=Found in Background 
c= !Found in Sequence	d= !Found in Background 

a = sequenceCounter		b= backgroundCounter
a+c=SeqsNumber			b+d=backgroundNumber
--> c= SeqsNumber - cumulated(sequenceCounter)
	d= backgroundNumber - cumulated(backgroundCounter)

One-sided Test:
	++a und --c ...
Two-sided Test:
	++a und --c
	--a und ++c ...



****************************/


void FisherExactTest(Seq &seq, 
					 Seq &back){

	
	double pValue=0;
	typedef std::map<String<Dna5>,unsigned int >::iterator MapIterator;
	MapIterator MapI=seq.seqCounter.begin();
	MapIterator MapIB=back.seqCounter.begin();
	
	for(;MapI!=seq.seqCounter.end();++MapI,++MapIB){
		

		modifyFET((*MapI).second,(*MapIB).second,(seq.SeqsNumber - (*MapI).second),(back.SeqsNumber - (*MapIB).second),pValue);
	
		seq.SortedPValue.insert(std::pair<double,String<Dna5> > (pValue, (*MapI).first));
		seq.SortedPValueReversed.insert(std::pair<String<Iupac>,double > ((*MapI).first,pValue));
	}


}

void FisherExactTest(std::map<String<Iupac>,unsigned int > &SequenceCounter,
					 std::map<String<Iupac>,unsigned int > &BackgroundCounter,
					 Seq &seq, 
					 Seq &back){

	
	

	double pValue=0;
	typedef std::map<String<Iupac>,unsigned int >::iterator MapIterator;
	MapIterator MapI=SequenceCounter.begin();
	MapIterator MapIB=BackgroundCounter.begin();
	
	for(;MapI!=SequenceCounter.end();++MapI,++MapIB){
		
		
		modifyFET((*MapI).second,(*MapIB).second,(seq.SeqsNumber - (*MapI).second),(back.SeqsNumber - (*MapIB).second),pValue);
	
		seq.generalizedSortedPValue.insert(std::pair<double,String<Iupac> > (pValue, (*MapI).first));
		seq.SortedPValueReversed.insert(std::pair< String<Iupac>,double> ((*MapI).first,pValue));
	}


}

/*
	Fisher-Exact-Test for generalized Kmere
*/
double FisherExactTest( Seq &seq,
						Seq &back,
						std::multimap<double,String<Iupac> > &GeneralizedSortedPValueTemp){
	
	
	
	if(seq.generalizedKmer.size()==0)	
		return 2;

	double pValue=0;
	typedef std::map<String<Iupac>,unsigned int >::iterator MapIterator;
	std::multimap<double,String<Iupac> >::iterator MapIterator2;
	MapIterator MapI = seq.generalizedKmer.begin();
	MapIterator MapIB= back.generalizedKmer.begin();
	

	for(;MapI!=seq.generalizedKmer.end();++MapI,++MapIB){
		
		if((*MapI).second > seq.SeqsNumber || (*MapIB).second > back.SeqsNumber){ std::cout<<(*MapI).first<<" "<<(*MapI).second<<" "<<(*MapIB).first<<" "<<(*MapIB).second<<std::endl;}
		modifyFET((*MapI).second,(*MapIB).second,seq.SeqsNumber - (*MapI).second,back.SeqsNumber - (*MapIB).second,pValue);
		

		GeneralizedSortedPValueTemp.insert(std::pair<double,String<Iupac> > (pValue, (*MapI).first));//not seq.generalizedSortedPValue, because this is the temp one
		seq.SortedPValueReversed.insert(std::pair<String<Iupac>,double > ((*MapI).first,pValue));
		}
	
	
	
	return GeneralizedSortedPValueTemp.begin()->first;
}




//void DebugFisherExact(unsigned int a,unsigned int b,unsigned int c,unsigned int d){
//	
//	double pValue=0;
//	if(c<a || d < b){
//		std::cerr<<"Cumulated Counter too large";
//		exit(1);
//	}
//	c=c-a;
//	d=d-b;
//	SEQAN_ASSERT_EQ((logFactorial(2+1+4+2)),(logFactorial(9)));
//	
//	std::cout<<"a: "<<a<<" b: "<<b<<" c: "<<c<<" d: "<<d<<std::endl;
//	pValue= logFactorial(a+b) + logFactorial(c+d) + logFactorial(a+c) + logFactorial(b+d) - (logFactorial(a+b+c+d) + logFactorial(a)+ logFactorial(b) + logFactorial(c) +logFactorial(d));
//	std::cout<<"log(pValue) "<<pValue<<std::endl;
//	pValue=logFactorial(a+b) + logFactorial(c+d) + logFactorial(a+c) + logFactorial(b+d) ;
//	std::cout<<"the dividend: "<<pValue<<std::endl;
//	pValue=- (logFactorial(2+1+4+2) + logFactorial(2)+ logFactorial(1) + logFactorial(4) +logFactorial(2));
//	std::cout<<"the divisor: "<<pValue<<std::endl;
//	pValue= (logFactorial(1));
//	std::cout<<"logFactorial(1): "<<pValue<<std::endl;
//	pValue= (logFactorial(0));
//	std::cout<<"logFactorial(0): "<<pValue<<std::endl;
//
//
//}

/****
		some Iupac-Mappings
		needed in the generalization step
****/
void MapIupac(IupacMaps &IMaps ){

	IMaps.IupacMap.get_allocator().allocate(16);
	IMaps.IupacMap[0]='U';
	IMaps.IupacMap[1]='T';
	IMaps.IupacMap[2]='A';
	IMaps.IupacMap[3]='W';
	IMaps.IupacMap[4]='C';
	IMaps.IupacMap[5]='Y';
	IMaps.IupacMap[6]='M';
	IMaps.IupacMap[7]='H';
	IMaps.IupacMap[8]='G';
	IMaps.IupacMap[9]='K';
	IMaps.IupacMap[10]='R';
	IMaps.IupacMap[11]='D';
	IMaps.IupacMap[12]='S';
	IMaps.IupacMap[13]='B';
	IMaps.IupacMap[14]='V';
	IMaps.IupacMap[15]='N';

	IMaps.IupacMapReversed.get_allocator().allocate(16);

	IMaps.IupacMapReversed['U']=0;
	IMaps.IupacMapReversed['T']=1;
	IMaps.IupacMapReversed['A']=2;
	IMaps.IupacMapReversed['W']=3;
	IMaps.IupacMapReversed['C']=4;
	IMaps.IupacMapReversed['Y']=5;
	IMaps.IupacMapReversed['M']=6;
	IMaps.IupacMapReversed['H']=7;
	IMaps.IupacMapReversed['G']=8;
	IMaps.IupacMapReversed['K']=9;
	IMaps.IupacMapReversed['R']=10;
	IMaps.IupacMapReversed['D']=11;
	IMaps.IupacMapReversed['S']=12;
	IMaps.IupacMapReversed['B']=13;
	IMaps.IupacMapReversed['V']=14;
	IMaps.IupacMapReversed['N']=15;

	

	IMaps.IupacMapReplace.get_allocator().allocate(14);
	IMaps.IupacMapReplace['R']="CT";//in Iupac-notation R=AG --> CT left
	IMaps.IupacMapReplace['Y']="AG";
	IMaps.IupacMapReplace['S']="AT";
	IMaps.IupacMapReplace['W']="CG";
	IMaps.IupacMapReplace['K']="AC";
	IMaps.IupacMapReplace['M']="GT";
	IMaps.IupacMapReplace['D']="C";
	IMaps.IupacMapReplace['H']="G";
	IMaps.IupacMapReplace['B']="A";
	IMaps.IupacMapReplace['V']="T";
	IMaps.IupacMapReplace['A']="CGT";
	IMaps.IupacMapReplace['C']="AGT";
	IMaps.IupacMapReplace['G']="ACT";
	IMaps.IupacMapReplace['T']="ACG";
	
	IMaps.IupacMapReplaceReversed.get_allocator().allocate(11);
	IMaps.IupacMapReplaceReversed['R']="AG";
	IMaps.IupacMapReplaceReversed['Y']="CT";
	IMaps.IupacMapReplaceReversed['S']="CG";
	IMaps.IupacMapReplaceReversed['W']="AT";
	IMaps.IupacMapReplaceReversed['K']="GT";
	IMaps.IupacMapReplaceReversed['M']="AC";
	IMaps.IupacMapReplaceReversed['D']="AGT";
	IMaps.IupacMapReplaceReversed['H']="ACT";
	IMaps.IupacMapReplaceReversed['B']="CGT";
	IMaps.IupacMapReplaceReversed['V']="ACG";
	IMaps.IupacMapReplaceReversed['N']="ACGT";
		
	IMaps.IupacMapInversed.get_allocator().allocate(11);
	IMaps.IupacMapInversed["AG"]='R';
	IMaps.IupacMapInversed["CT"]='Y';
	IMaps.IupacMapInversed["CG"]='S';
	IMaps.IupacMapInversed["AT"]='W';
	IMaps.IupacMapInversed["GT"]='K';
	IMaps.IupacMapInversed["AC"]='M';
	IMaps.IupacMapInversed["AGT"]='D';
	IMaps.IupacMapInversed["ACT"]='H';
	IMaps.IupacMapInversed["CGT"]='B';
	IMaps.IupacMapInversed["ACG"]='V';
	IMaps.IupacMapInversed["ACGT"]='N';


}
/*	- initiate by repeatedly picking the top motifs from SortedPValue
	- calls GeneralizeKmer and the EstimateFunction
*/
void InitGeneralization(IupacMaps &IMaps,
						Seq &seq,
						Seq &back){
	unsigned int i=0;
	unsigned int limit;
	std::multimap<double,String<Dna5> >::iterator MapIterator;
	std::multimap<double,String<Iupac> >::iterator MapIteratorT;	
	std::multimap<double,String<Iupac> > generalizedSortedPValueTemp;
	std::map<String<Iupac>,unsigned int> generalizedKmerTemp;
	std::map<String<Iupac>,unsigned int> generalizedKmerBackgroundTemp;
	
	if(seq.SortedPValue.size()>seq.seed)	limit=seq.seed;//standard-seed = 100
	else if(seq.SortedPValue.size()==0) return;
	else	limit = seq.SortedPValue.size();
	
	for(MapIterator=seq.SortedPValue.begin();i<limit;++MapIterator,++i){//iterate over Top100(TopSeedNumber)
		GeneralizeKmer((*MapIterator).second,IMaps,seq,back);
	}
	/*
		- only do the next function call, if in the last at least one pValue<treshold	
		- call GeneralizeKmer in loop
	*/
	
	double topPValue = FisherExactTest(seq,back,generalizedSortedPValueTemp);// lowest pValue from the first generalization
	double topPValueOld =seq.SortedPValue.begin()->first;//lowest pValue before generalization
	
	while(topPValue<0.05 && topPValue<topPValueOld){//only start a new round, if the top PValue is an improvement of the old one
		
		
		/*
			while wird das erste mal mit generalizedKmer aufgerufen und dem temporären mapping der pValues
			das temporäre mapping wird in das richtige mapping gemerged und gecleant, damit geschaut werden kann, ob bei den neuen pValues ein wert
			über dem treshold ist --> falls nicht bricht die while ab
			falls doch wird generalizedKmer kopiert und gecleant aus dem gleichen grund, damit nur die neuen generalisierten Kmere untersucht werden
			--> das Temp hier, um über alle alten zu gehen und um diese weiter zu generalisieren

		*/

		seq.generalizedSortedPValue.insert(generalizedSortedPValueTemp.begin(),generalizedSortedPValueTemp.end());

	
		generalizedKmerTemp.clear();
		generalizedKmerBackgroundTemp.clear();
		generalizedKmerTemp=seq.generalizedKmer;
		generalizedKmerBackgroundTemp=back.generalizedKmer;
		back.generalizedKmer.clear();
		seq.generalizedKmer.clear();
		if(generalizedSortedPValueTemp.size()>seq.seed)	limit=seq.seed;
		else if(generalizedSortedPValueTemp.size()==0) return;
		else	limit = generalizedSortedPValueTemp.size();
	
		i=0;//only Top100 (seed number)
		for(MapIteratorT=generalizedSortedPValueTemp.begin();i<limit;++MapIteratorT,++i){//iterate over Top100 (seed number)
		
			//Temp ums zu finden, aber das normale auch übergeben, zum neu befüllen
			
			GeneralizeKmer((*MapIteratorT).second,generalizedKmerTemp,generalizedKmerBackgroundTemp,IMaps,seq,back);
			
			                                                                        
	
		}
	generalizedSortedPValueTemp.clear();
	topPValueOld =topPValue;
	topPValue = FisherExactTest(seq,back,generalizedSortedPValueTemp);
		
	};
	
	
}
/*	- gets a Kmer and replaces each position successively with each possible ambiguity Code(Iupac)
	- only one wildcard per String at one time
	- the unsigned int corresponds to the estimated counter
	- String<Dna5> for initialization
 */
void GeneralizeKmer(String<Dna5> Kmer,
					IupacMaps &IMaps,
					Seq &seq,
					Seq &back){
	
	String<Iupac> temp;//temporary String --> generalizedKmer[temp]
	String<Iupac> temp2;//Kmer with replaced position--> relevant for estimateCounter
	Iterator<String<Iupac> >::Type tempIt;//Iterator over temp --> same length as Kmer
	String<Dna5> replace = "ACGT";//replace the current position with each possible ambiguity code --> A,C,G or T
	Iterator<String<Dna5> >::Type replaceIt;
	unsigned int counter =0;
	char tempChar;
	//replaceIt = begin(replace);	
	temp = Kmer;
	temp2=Kmer;
	tempIt = begin(temp);
	
	
	for(;tempIt!=end(temp);++tempIt){//loop over each position in kmer
		replaceIt = begin(replace);	
		for(;replaceIt!=end(replace);++replaceIt){// loop over ACGT
			temp = Kmer;// reset temp
			
			if(*tempIt == *replaceIt) continue; // there is no Iupac for "AA", in return "AG" = "R"
			
			tempChar =*tempIt;//stores the current char because of temp2
			*tempIt=*replaceIt;
			temp2=temp;//temp2 stores the char for the new wildcard
			*tempIt=tempChar;//temp ist the temp from before -->gets expanded with a new wildcard in the next step
			
			*tempIt =IMaps.IupacMap[IMaps.IupacMapReversed[*tempIt] + IMaps.IupacMapReversed[*replaceIt]];//compute Iupac-letter--> A + G = R and replace the current location in temp
			
			if(seq.generalizedKmer.find(temp)!=seq.generalizedKmer.end()) continue;// if Kmer is in the Map -->nothing to do
			//call estimateCounter with Kmer and temp2 --> Kmer=AAA  temp2=TAA temp=WAA
			estimateCounter(seq,Kmer,temp2,counter);
			
			seq.generalizedKmer[temp]=counter;//temp = new motif
			estimateCounter(back,Kmer,temp2,counter);
			back.generalizedKmer[temp]=counter;
		}
	}


	

}

/*	- the same as above except that each String has already a wildcard
*/
void GeneralizeKmer(String<Iupac> Kmer,
					std::map<String<Iupac>,unsigned int> &generalizedKmerTemp,
					std::map<String<Iupac>,unsigned int> &generalizedKmerBackgroundTemp,
					IupacMaps &IMaps,
					Seq &seq, 
					Seq &back){
	unsigned counter =0;
	char tempChar;
	String<Iupac> temp;//temporary String --> generalizedKmer[temp]
	String<Iupac> temp2;//Kmer with replaced position--> relevant for estimateCounter
	Iterator<String<Iupac> >::Type tempIt;//Iterator over temp --> same length as Kmer
	String<Iupac> replace;
	Iterator<String<Iupac> >::Type replaceIt;
	
	
	temp = Kmer;
	tempIt = begin(temp);
	for(;tempIt!=end(temp);++tempIt){//loop over each position in kmer
				
		if(*tempIt =='N')continue;//nothing to do
		replace=IMaps.IupacMapReplace[*tempIt];
		replaceIt = begin(replace);	
		for(;replaceIt!=end(replace);++replaceIt){// loop over the replacement-chars, W=AT --> replace=CG
			temp = Kmer;// reset temp
			tempChar =*tempIt;//stores the current char because of temp2
			*tempIt=*replaceIt;//replace the current char for temp2
			temp2=temp;
			*tempIt=tempChar;

			*tempIt =IMaps.IupacMap[IMaps.IupacMapReversed[*tempIt] + IMaps.IupacMapReversed[*replaceIt]];
			
			if(seq.SortedPValueReversed[Kmer] >= 0.05 || seq.SortedPValueReversed[temp2] >= 0.05) continue;//only if Kmer and temp2 are significant estimate the counter
			if(generalizedKmerTemp.find(temp)!=generalizedKmerTemp.end()) continue;
			estimateCounter(seq,generalizedKmerTemp,Kmer,temp2,counter);
			seq.generalizedKmer[temp]=counter;
			estimateCounter(back,generalizedKmerBackgroundTemp,Kmer,temp2,counter);
			back.generalizedKmer[temp]=counter;
		}

		}

}


//void loopOverRecplacement(String<Iupac> &temp,String<Iupac> &temp2,String<Iupac> Kmer,Iterator<String<Iupac> >::Type tempIt,unsigned int counter){
//	char tempChar;
//	String<Iupac> replace;
//	Iterator<String<Iupac> >::Type replaceIt;
//
//	replace=IupacMapReplace[*tempIt];
//	replaceIt = begin(replace);	
//
//	for(;replaceIt!=end(replace);++replaceIt){// loop over the replacement-chars
//			temp = Kmer;// reset temp
//			tempChar =*tempIt;//stores the current char because of temp2
//			*tempIt=*replaceIt;//replace the current char for temp2
//			temp2=temp;
//			*tempIt=tempChar;
//
//			*tempIt =IupacMap[IupacMapReversed[*tempIt] + IupacMapReversed[*replaceIt]];
//			//only if Kmer and temp2 are significant estimate the counter
//			if(generalizedKmerTemp.find(temp)!=generalizedKmerTemp.end()) continue;
//			estimateCounter(SequenceCounter,generalizedKmerTemp,Kmer,temp2,counter,SeqsNumber);
//			generalizedKmer[temp]=counter;
//			estimateCounter(BackgroundCounter,generalizedKmerBackgroundTemp,Kmer,temp2,counter,BackgroundNumber);
//			generalizedKmerBackground[temp]=counter;
//			//std::cout<<temp<<" "<<counter<<std::endl;
//		}
//
//}




/*
	- estimates the Counter for the initial wildcard-Kmere
	- SequencesCounter and BackgroundCounter
*/
void estimateCounter(Seq &seq,
					 String<Iupac> temp,
					 String<Iupac> temp2,
					 unsigned int &counter){
	unsigned int RE1=seq.seqCounter.find(temp)->second;//old motif -->from the last step
	
	if(seq.seqCounter.find(temp2)!=seq.seqCounter.end()){//if temp2=known motif --> counter available
		counter= RE1+ seq.seqCounter.find(temp2)->second - (RE1*seq.seqCounter.find(temp2)->second)/seq.SeqsNumber;
		if(counter>seq.SeqsNumber){ std::cerr<<"Error, counter to big in if"<<std::endl;
		std::exit(1);
		}
	}
	else{
		counter=RE1;//RE2=0
		if(counter>seq.SeqsNumber){ std::cerr<<"Error, counter to big in else"<<std::endl;
		std::exit(1);
		}
	}
	
}
/*
	- estimated the Counter for the next Kmer
*/
void estimateCounter(Seq &seq,
					 std::map<String<Iupac>,unsigned int> &generalizedKmer,
					 String<Iupac> temp,
					 String<Iupac> temp2,
					 unsigned int &counter){


	if(generalizedKmer.find(temp)== generalizedKmer.end()){
		std::cerr<<"Error, could not find "<<temp<<" in generalizedKmer"<<std::endl;
		std::exit(1);
	}

	
	
	unsigned int RE1=(*generalizedKmer.find(temp)).second;//the new seed RE is a Kmer with wildcard
	//temp2 may be in generalizedKmer or in SequenceCounter
	if(seq.seqCounter.find(temp2)!=seq.seqCounter.end()){// if temp2 is in SequenceCounter do the same as above --> has no wildcard
		counter= RE1+ seq.seqCounter.find(temp2)->second- (RE1*seq.seqCounter.find(temp2)->second)/seq.SeqsNumber;
		if(counter>seq.SeqsNumber){ std::cerr<<"Error, counter to big in if"<<std::endl;
		std::exit(1);
		}
	}
		
	else if(generalizedKmer.find(temp2)!=generalizedKmer.end()){//if temp2 has a wildcard and is found in generalizedKmer
		counter= RE1+ generalizedKmer.find(temp2)->second - (RE1*generalizedKmer.find(temp2)->second)/seq.SeqsNumber;
		if(counter>seq.SeqsNumber){ std::cerr<<"Error, counter to big in elseif"<<std::endl;
		std::exit(1);
		}
		
	}
	else{//if temp2 is not found
		counter= RE1;//RE2=0
		if(counter>seq.SeqsNumber){ std::cerr<<"Error, counter to big in else"<<std::endl;
		std::exit(1);
		}
	}
}

/****
		count the generalized Kmer exact
****/
void exactGeneralizeCount(  std::map<String<Iupac>,unsigned int > &seqCounter,
							std::map<String<Iupac>,unsigned int > &backCounter,
							Finder<Index<StringSet<String<Dna5> > > > &finder,
							Finder<Index<StringSet<String<Dna5> > > > &finderB,
							Seq &seq,
							Seq &back,
							IupacMaps &IMap){

	std::multimap<double,String<Iupac> >::iterator generalizedSortedPValueIt;

	generalizedSortedPValueIt = seq.generalizedSortedPValue.begin();

	for(unsigned int i=0;i<seq.seed && generalizedSortedPValueIt!=seq.generalizedSortedPValue.end() ;++i,++generalizedSortedPValueIt){
		if(seqCounter.find((*generalizedSortedPValueIt).second)!=seqCounter.end()) continue;
		
		CountKmer(seqCounter,finder,(*generalizedSortedPValueIt).second,seq,IMap);
		CountKmer(backCounter,finderB,(*generalizedSortedPValueIt).second,back,IMap);
		
	}
	
	seq.generalizedSortedPValue.clear();
	
	FisherExactTest(seqCounter,backCounter,seq,back);//computes the pValue of each Motif due to the counter
	seqCounter.clear();
	backCounter.clear();



}


/****
		computes the alginment of two pwms with a slightly modified Smith-Waterman-Algorithm
****/
String<double> AlignPWMs(	std::map<unsigned int,std::map<Iupac,double> > &freqMatrix1,
							std::map<unsigned int,std::map<Iupac,double> > &freqMatrix2){

  int freqL1 = length(freqMatrix1);                     
  int freqL2 = length(freqMatrix2);

  std::vector<std::vector<double> > M(freqL1+1,std::vector<double>(freqL2+1));     
  for(unsigned i=0;i<=freqL1;++i){
    for(unsigned j=0;j<=freqL2;++j){
      M[i][j]=0;
    }
  } 
 
  
	//no traceback needed-->no gaps allowed, except for the opening-gaps
	String<double> Mmax;
	appendValue(Mmax,100);
	appendValue(Mmax,0);
	appendValue(Mmax,0);

	
	for(unsigned i=1;i<=freqL1;++i){
		for(unsigned j=1;j<=freqL2;++j){

			M[i][j]=M[i-1][j-1]+ComparePWM(freqMatrix1[i-1],freqMatrix2[j-1],Entropy());//the more different, the bigger is the value
			
			
				
			if((i==freqL1 && j>freqL2*0.5  && M[i][j]/j<Mmax[0])){//normalize to the length of the overlap--> /j --> Overlap must be greater than 0.5
				Mmax[0]=M[i][j]/j;
				Mmax[1]=i;
				Mmax[2]=j;
			}
			else if((j==freqL2 &&i>freqL1*0.5 && M[i][j]/i<Mmax[0])){
				Mmax[0]=M[i][j]/i;
				Mmax[1]=i;
				Mmax[2]=j;
				
			}
			
			
		}
		
		
	}
	
	return Mmax;
}



#endif  // #ifndef SANDBOX_MEYERCLP_APPS_DREME_H_
