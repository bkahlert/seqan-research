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

using namespace seqan;


struct Seq
{
	StringSet<CharString> ids;
	unsigned c; // Anzahl der Motive, provisorisch
	StringSet<String<Dna5> > seqs;//
	Index< StringSet<String<Dna5> > > SArray;
	unsigned seed;
	unsigned int SeqsNumber;
	
	std::map<String<Dna5>,unsigned int > seqCounter;//maps the Sequence-Kmere to a Counter for the Sequence
	std::map<String<Iupac>,unsigned int> generalizedKmer;

	std::map<String<Iupac>,double > SortedPValueReversed;
	
	std::multimap<double,String<Dna5> > SortedPValue;
	std::multimap<double,String<Iupac> > generalizedSortedPValue;
	std::map<unsigned int,std::map<Iupac,double> > freqMatrix;
	String< std::map<unsigned int,std::map<Iupac,double> > > allPWMs;
	std::map<unsigned int,std::map<Iupac,double> > weightMatrix;
	std::map<unsigned int, double > seqLogoMatrix;
	std::map<unsigned int, double > columnEntropy;
	
	std::map<unsigned int,std::map<Iupac,double> > InformationContentMatrix;
	FrequencyDistribution<Dna5> frequencies;

	/**
		In BuildFrequencyMatrix werden die gefundenen Intervalle des Top-Motivs an Intervalls gehängt
		Am Ende der while wird dann ein IntervallTree für jede Sequenz erzeugt, sodass im nächsten Schritt direkt geschaut werden kann ob ein Motiv mit einem schon 
		gefundenen überlappt
	**/
	typedef IntervalAndCargo<unsigned, unsigned> TInterval;
	String<String<TInterval> > intervals;//String<TInterval> enthält alle Intervalle einer Sequenz. String<String<..>> enthält alle Sequenzen
	
	String<IntervalTree<unsigned> > intervalTrees; // Tree für schnellere Suche erstellen --> String von Trees, da mehrere Sequenzen
	String<unsigned> results;
	double pValue;
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


struct Euklidisch_;
typedef Tag<Euklidisch_> Euklidisch;

struct Pearson_;
typedef Tag<Pearson_> Pearson;

struct Mahalanobis_;
typedef Tag<Mahalanobis_> Mahalanobis;

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
double ComparePWM(Seq &seq,std::map<unsigned int,std::map<Iupac,double> > &freqMatrix1,std::map<unsigned int,std::map<Iupac,double> > &freqMatrix2, Entropy const & tag);
double ComparePWM(Seq &seq,std::map<unsigned int,std::map<Iupac,double> > &freqMatrix1,std::map<unsigned int,std::map<Iupac,double> > &freqMatrix2, Euklidisch const & tag);
void PWMClustering(String< std::map<unsigned int,std::map<Iupac,double> > > &allPWMs, String<Cluster> &cluster);
void UpdateDistantMatrix(int n, int x, int y, std::vector<std::vector<String<double > > > &compare, CompleteLinkage const & tag);
void UpdateDistantMatrix(int n, int x, int y, std::vector<std::vector<String<double > > > &compare,std::vector<unsigned> &weights, AverageLinkage const & tag);
void minDifferenceInMatrix(unsigned n,String<double> &minDifference,std::vector<std::vector<String<double > > > compare);
void replaceKmer(Seq &seq,unsigned int stringNumber, unsigned int begin, unsigned int end);
void saveData(Seq &seq,std::ofstream &PWM,unsigned c);
void computesDistantMatrix(String< std::map<unsigned int,std::map<Iupac,double> > > &allPWMs,std::vector<std::vector<String<double > > > &compare, unsigned allPWMsLength);
String<double> AlignPWMs(std::map<unsigned int,std::map<Iupac,double> > &freqMatrix1,std::map<unsigned int,std::map<Iupac,double> > &freqMatrix2);
double computeDr(Cluster &cluster,std::vector<std::vector<String<double> > > copyCompare);
double computeWk(int n,String<Cluster> &cluster);
void computeReferenceData(String< std::map<unsigned int,std::map<Iupac,double> > > &Reference,unsigned allPWMsLength,unsigned PWMLength );
void computeGapStat(String< std::map<unsigned int,std::map<Iupac,double> > > &allPWMs, String<Cluster> &observedCluster);
void compute_l_quer(String<Cluster> &cluster, String<double> &l_quer_for_k,String<String<double> > &allWk,unsigned i,int B,	unsigned allPWMsLength);
	void computeGap(String<double> &l_quer_for_k,String<Cluster> &observedCluster,String<double> &Gap,unsigned allPWMsLength);
void compute_sdk_and_sk(String<double> &sk,String<double> &l_quer_for_k,String<String<double> > &allWk,unsigned allPWMsLength,int B);


void readFastA( struct Seq &seq, 
				CharString fname){

	//########### einlesen in ids und seqs
	//Variable um Sequenz und ID zu speichern
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

	typedef Index<StringSet<String<Dna5> > > TMyIndex;//Dna5Q
	Iterator<TMyIndex, BottomUp<> >::Type myIterator(seq.SArray);
	for(;!atEnd(myIterator);++myIterator){
		std::cout<<representative(myIterator)<<std::endl;
		
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


//iniate Search in Fore- and Background
void initExactKmer( Seq &seq,
					Seq &back,
					unsigned int kmer_len,
					unsigned int kmer_len_end){

	Finder<Index<StringSet<String<Dna5> > > > finder(seq.SArray);
	Finder<Index<StringSet<String<Dna5> > > > finderB(back.SArray);//finder background
	typedef Index< StringSet<String<Dna5> > > TMyIndex;
	//kmer_len= 3;//minimal kmer-length
	//kmer_len_end=8;//maximal length
	if(kmer_len<1) kmer_len=3;
	if(kmer_len_end<kmer_len) kmer_len_end=kmer_len+1;
	//std::cout<<"kmer: "<<kmer_len<<std::endl;
	//std::cout<<"end: "<<kmer_len_end<<std::endl;

	typedef Iterator<StringSet<String<Dna5> > >::Type TStringSetIterator;
	unsigned int slen=0;
	String<Dna5> Kmer;//current Kmer
	

	std::cout<<std::endl<<std::endl;
	for(;kmer_len<=kmer_len_end;++kmer_len){//loop over all possible Kmer-length -->3-8
		for(TStringSetIterator it=begin(seq.seqs);it!=end(seq.seqs);++it){//loop over all sequences
			slen=length(value(it));//length of the current seq
			if(slen==0) continue;
			else if(slen<kmer_len) continue;
			//std::cout<<"Sequence: "<<value(it)<<std::endl;
		
		
			for(unsigned int i=0;i<slen-kmer_len+1;++i){//loop over all Kmere in the current sequence
			
				Kmer=infix(value(it),i,i+kmer_len);//optimieren
				//std::cout<<Kmer<<" ";
				//if(Kmer[kmer_len-1]=='N'){//'AAAN'AAAA ---> AAA'NAAA'A --> after continue AAAN'AAAA'
				//	i+=kmer_len;
				//	continue;
				//}
				
				if(seq.seqCounter.find(Kmer)!=seq.seqCounter.end()) continue;// if Kmer is in the Map -->nothing to do
	//			//Pattern<String<Dna5> > pattern(Kmer);
			//	std::cout<<"count";
				CountKmer(back,finderB,Kmer);
				CountKmer(seq,finder,Kmer);
				
			}
		}

	}

		
	//}
}

//gets the current Kmer and searches it in the index
//max(cumulated sum of the counter)=SeqsNumber --> counts number of sequences containing the motif
void CountKmer( Seq &seq, 
				Finder<Index<StringSet<String<Dna5> > > > &finder, 
				String<Dna5> &Kmer){
	
			std::vector<int> CounterV(seq.SeqsNumber+1,0);//counter for storing 1 or 0 for each Seq + the cumulated sum of the counter in the last field
			//std::cout<<"vor while  ";
			clear(finder);
			while(find(finder,Kmer)){//search the current Kmer in all sequences
				//std::cout<<'[' <<beginPosition(finder)<<','<<endPosition(finder)<<")\t"<<infix(finder)<<std::endl;//Debug
			
				if(seq.c>1){//nur im foreground
					clear(seq.results);
					findIntervals(seq.intervalTrees[beginPosition(finder).i1], beginPosition(finder).i2, endPosition(finder).i2, seq.results);
				}// wenn results>0, dann überlappt das Kmer --> nicht aufzählen
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
		//if not end call fkt. with temp 
		// if end call find --> &counter
		tempIttemp=tempIt;//der rekursive aufruf mit diesem, da die schleife mit tempIt weitergehen soll
		if(tempIt+1!=end(temp)){
			loopOverKmer(seq,temp,temp,++tempIttemp,finder,counter,CounterV,IMap);

		}
		
		if((*tempIttemp == 'A' || *tempIttemp == 'C' || *tempIttemp == 'G' || *tempIttemp == 'T' || tempIttemp == end(temp))){
					
			
			FindKmer(seq,temp,finder,counter,CounterV);
		}
				
		*tempIt=resetTemp;			
					
	}
				
			//}



}

void FindKmer(  Seq &seq,
				String<Iupac> &temp,
				Finder<Index<StringSet<String<Dna5> > > > &finder,
				unsigned int &counter,
				std::vector<int> &CounterV){



					clear(finder);
					while(find(finder,temp)){//search the current Kmer in all sequences
						//std::cout<<'[' <<beginPosition(finder)<<','<<endPosition(finder)<<")\t"<<infix(finder)<<std::endl;//Debug
						
						if(seq.c>1){//nur im foreground
							clear(seq.results);
							findIntervals(seq.intervalTrees[beginPosition(finder).i1], beginPosition(finder).i2, endPosition(finder).i2, seq.results);
						}
						if(CounterV[beginPosition(finder).i1] == 0 && length(seq.results)==0){//count number of sequences containing the motif, not the occurrences to avoid problems with self-overlapping
							//ansonsten m√ºsste das array noch einmal durch gegangen werden und an jeder stellt !=0 ++
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
			
			//std::cout<<"Kmer "<<Kmer<<std::endl;
			loopOverKmer(seq,temp,Kmer,tempIt,finder,counter,CounterV,IMap);

			
			
			seqCounter[Kmer]=counter;
			//std::cout<<Kmer<<" "<<seqCounter[Kmer]<<"   ";
			//std::cout<<temp<<" "<<CounterV[SeqsNumber]<<std::endl;
			CounterV.clear();
			
			/*
				Der clear hier bewirkt, dass bei ASG entweder ACG oder AGG vorkommen darf,
				falls beide vorkommen(in einer Sequenz) wird der counter trotzdem nur um 1 erh√∂ht
				-->Counter f√ºr AGG ist falsch, ASG stimmt jedoch wieder. Counter[AGG] --> irrelevant
				-->AGG wird jedes mal wenns ben√∂tigt wird neu berechnet-->optimieren
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
	//std::cout<<temp<<" vor-while"<<std::endl;
	while(find(finder,temp)){//search the current Kmer in all sequences
		if(seq.c==1){//nur im ersten Schritt in Intervals speichern, danach kann direkt an den Tree angefügt werden
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
	for( unsigned int k =0;k< length(temp);++k){//berechnet Frequenz der Nukleotide, wobei jedes Motiv wieder nur einmal pro Sequenz zählt!
		seq.freqMatrix[k][temp[k]]+=CounterV[seq.SeqsNumber]; //GCAGCA --> counter der einzelnen wird um die gleiche anzahl hochgezählt, GCAGTA --> usw. Nicht-Wildcards haben W'keit 1
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
	//std::cout<<temp<<" vor-while"<<std::endl;
	while(find(finder,temp)){//search the current Kmer in all sequences
		if(seq.c==1){//nur im ersten Schritt in Intervals speichern, danach kann direkt an den Tree angefügt werden
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
	for( unsigned int k =0;k< length(temp);++k){//berechnet Frequenz der Nukleotide, wobei jedes Motiv wieder nur einmal pro Sequenz zählt!
		seq.freqMatrix[k][temp[k]]+=CounterV[seq.SeqsNumber]; //GCAGCA --> counter der einzelnen wird um die gleiche anzahl hochgezählt, GCAGTA --> usw. Nicht-Wildcards haben W'keit 1
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
			
			char resetTemp;//bei mehr als einer wildcard, müssen die weiter hinten liegenden nach abarbeitung resetet werden, ansonsten werden diese im nächsten schritt übergangen

			if(tempIt==end(temp)) return;//&&(tempIt+1!=end(temp))
			/*freq[*tempIt]=1;
			freqMatrix[pos]=freq;
			freqMatrix[pos]['A']=1;
			freq.clear();*/
			//std::cout<<" "<<*tempIt<<" ";
			if((*tempIt == 'A' || *tempIt == 'C' ||*tempIt == 'G' ||*tempIt == 'T')){
				loopOverTopKmer(seq,temp,temp,++tempIt,finder,counter,CounterV,IMaps);//only replace the position with a wildcard
				return;//nach diesem schritt immer return, sonst gelangt man in eine loop
			}
				
			
				replace=IMaps.IupacMapReplaceReversed[*tempIt];
				replaceIt = begin(replace);	
				for(;replaceIt!=end(replace);++replaceIt){
					//std::cout<<" "<<temp<<" "<<Kmer<<" "<<*replaceIt<<std::endl;
					temp = Kmer;// reset temp
					resetTemp = IMaps.IupacMapInversed[replace]; //falls Y ersetzt wird, ist replace CT --> also resetTemp wieder Y
					//std::cout<<" resetTemp "<<resetTemp<<std::endl;
					*tempIt = *replaceIt;
					//std::cout<<" re "<<temp<<" ";
					tempIttemp=tempIt;//der rekursive aufruf mit diesem, da die schleife mit tempIt weitergehen soll
				//	std::cout<<" "<<temp<<" "<<Kmer<<" "<<*replaceIt<<std::endl;
					
					if(tempIt+1!=end(temp)){
						//std::cout<<"vor if "<<temp<<std::endl;
						loopOverTopKmer(seq,temp,temp,++tempIttemp,finder,counter,CounterV,IMaps);
						//std::cout<<*tempIttemp<<" tempittemp "<<std::endl;
					}
					//zu häufig aufgerufen
					if((*tempIttemp == 'A' || *tempIttemp == 'C' || *tempIttemp == 'G' || *tempIttemp == 'T' || tempIttemp == end(temp))){//falls nicht, dann kann dies übersprungen werden
						/**
						posTemp gibt die Stelle an, an der man sich grad befindet 
						--> freqMatrix enthält für jede Position die W'keit für ACGT
						--> counter ist die Anzahl wie oft das Wildcard Motiv insgesamt gefunden wurde
						--> CounterV an letzter Stelle die Anzahl für das jeweilige nicht-Wildcard Motiv
						--> bei GSAKYA z.B. als Motiv wird jedes Motiv bei 'S' vier mal gesucht(durch die anderen 2 Wildcards)
						--> CounterV für 1 bis posTemp aufaddieren --> in freqMatrix und zwar für die jeweiligen *tempIt-chars
						--> am Ende alle durch counter teilen --> aufpassen, für jeweilige pos gibts verschiedene counter
						--> FindKmer wird nur mit ganzen aufgerufen, also alle addieren, dann ist der counter auch gleich?
						**/
						std::vector<int> CounterV(seq.SeqsNumber+1,0);
						FindTopKmer(seq,temp,finder,counter,CounterV);
						
					}
					*tempIt=resetTemp;
					//if ende von replaceIt, dann begin(replace) in temp speichern und mitübergeben--> referenz?
				}


}

/***
	Computes PWM 
***/
void BuildFrequencyMatrix(  Finder<Index<StringSet<String<Dna5> > > > &finder,
							String<Iupac> &Kmer,
							Seq &seq, 
							IupacMaps &IMaps){
			//std::cout<<Kmer<<std::endl;
			//freqMatrix -->unsigned int = position in Kmer, position 1 in map = prob. for A, pos. 2 = prob. for C...
			String<Iupac> temp;	

			Iterator<String<Iupac> >::Type tempIt;

			temp=Kmer;
			tempIt = begin(temp);
			unsigned int counter=0;
			
			std::vector<int> CounterV(seq.SeqsNumber+1,0);
		

			loopOverTopKmer(seq,temp,Kmer,tempIt,finder,counter,CounterV,IMaps);
			CounterV.clear();
			//loopOver funktionier, aber jetzt wird der counter nicht mehr richtig berechnet --> fixen  + andere loop anpassen
			//Durch die wildcards mehrere Vorkommen pro Sequence möglich:
			//seqCounter[temp]=CounterV[seq.SeqsNumber];
					
			//std::cout<<temp<<" "<<*replaceIt<<" ";
			
					
					
			
				
			if(counter>0){//normalisieren des Counters
				
				
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
		Version für nicht generalisierte
**/
void BuildFrequencyMatrix(  Finder<Index<StringSet<String<Dna5> > > > &finder,
							String<Dna5> &Kmer,
							Seq &seq, 
							IupacMaps &IMaps){
			//std::cout<<Kmer<<std::endl;
			//freqMatrix -->unsigned int = position in Kmer, position 1 in map = prob. for A, pos. 2 = prob. for C...
			String<Iupac> temp;	
			
			Iterator<String<Iupac> >::Type tempIt;
			
			temp=Kmer;
			tempIt = begin(temp);
			unsigned int counter=0;
			
			std::vector<int> CounterV(seq.SeqsNumber+1,0);
			FindTopKmer(seq,Kmer,finder,counter,CounterV);
			CounterV.clear();
			
			
				
			if(counter>0){//normalisieren des Counters
				
				
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

void BuildWeightMatrix(Seq &seq){

	
	
	for(unsigned int j=0;j<length(seq.freqMatrix);++j){
		seq.weightMatrix[j]['A'] = log(seq.freqMatrix[j]['A']/seq.frequencies[0]);
		seq.weightMatrix[j]['C'] = log(seq.freqMatrix[j]['C']/seq.frequencies[1]);
		seq.weightMatrix[j]['G'] = log(seq.freqMatrix[j]['G']/seq.frequencies[2]);
		seq.weightMatrix[j]['T'] = log(seq.freqMatrix[j]['T']/seq.frequencies[3]);
	}




}

void saveData(Seq &seq,std::ofstream &PWM,unsigned c){

	String<char> pwm;
	char *buffer= new char[33];
	append(pwm,"C:/Users/David/Desktop/PWM");
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
	PWM.close();
	delete[] buffer;

}

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

double ComparePWM(std::map<Iupac,double>  &freqMatrix1,std::map<Iupac,double>  &freqMatrix2, Entropy const & tag){


	
	/****
			Für die übergebenen Spalten der Matrizen wird die Entropy berechnet
			Eintrag ist = 0, wenn die Werte identisch sind. Je größer die Zahl, desto unterschiedlicher die Werte
			Was wenn Wert = Hintergrundverteilung? --> Eintrag wäre = 0, obwohl das nichts mit dem Motiv zu tun hat

			Für jede Spalte berechnen und durch Spaltenanzahl teilen --> im Nachhinein
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

double ComparePWM(std::map<Iupac,double>  &freqMatrix1,std::map<Iupac,double>  &freqMatrix2, Euklidisch const & tag){

	double columnEntropy = 0;

	columnEntropy  = (freqMatrix1['A'] - freqMatrix2['A'])*(freqMatrix1['A'] - freqMatrix2['A']);
	columnEntropy += (freqMatrix1['C'] - freqMatrix2['C'])*(freqMatrix1['C'] - freqMatrix2['C']);
	columnEntropy += (freqMatrix1['G'] - freqMatrix2['G'])*(freqMatrix1['G'] - freqMatrix2['G']);
	columnEntropy += (freqMatrix1['T'] - freqMatrix2['T'])*(freqMatrix1['T'] - freqMatrix2['T']);
	columnEntropy  = sqrt(columnEntropy);

	return columnEntropy;
	
}
/****
		Speichert die jeweiligen Mittelwerte der 2 übergebenen Matrizen 
****/
void BuildMeanOf2PWMs(Seq &seq,std::map<unsigned int,std::map<Iupac,double> > &freqMatrix1,std::map<unsigned int,std::map<Iupac,double> > &freqMatrix2){

	for(unsigned int j=0;j<length(freqMatrix1);++j){
		freqMatrix1[j]['A'] = (freqMatrix1[j]['A']+freqMatrix2[j]['A'])/2;
		freqMatrix1[j]['C'] = (freqMatrix1[j]['C']+freqMatrix2[j]['C'])/2;
		freqMatrix1[j]['G'] = (freqMatrix1[j]['G']+freqMatrix2[j]['G'])/2;
		freqMatrix1[j]['T'] = (freqMatrix1[j]['T']+freqMatrix2[j]['T'])/2;
	}
	

}

void UpdateDistantMatrix(int n, int x, int y, std::vector<std::vector<String<double> > > &compare, CompleteLinkage const & tag ){

	/****
			Updates the distance with Complete Linkage
			Replaces PWM x with the new one
			
	****/
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




void UpdateDistantMatrix(int n, int x, int y, std::vector<std::vector<String<double> > > &compare,std::vector<unsigned> &weights, AverageLinkage const & tag ){


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


void minDifferenceInMatrix(unsigned n,String<double> &minDifference,std::vector<std::vector<String<double> > > compare){

	minDifference[0]=0;
	for(unsigned i=0;i<n-1;++i){

		for(unsigned j=i+1;j<n;++j){
			
			if(minDifference[0]==0 || minDifference[0]>compare[j][i][0]){
				minDifference[0]=compare[j][i][0];
				minDifference[1]=i;
				minDifference[2]=j;
				//std::cout<<i<<" "<<j<<" "<<compare[j][i][0];
			}

		}

	}

	//std::cout<<std::endl<<minDifference[0]<<" "<<minDifference[1]<<" "<<minDifference[2]<<std::endl<<std::endl;
}


void computesDistantMatrix(String< std::map<unsigned int,std::map<Iupac,double> > > &allPWMs,std::vector<std::vector<String<double> > > &compare, unsigned allPWMsLength){

	unsigned j;
	unsigned i;
	for(i=0;i<allPWMsLength-1;++i){//computes distances with ComparePWM and saves it in the matrix compare

		for( j=i+1;j<allPWMsLength;++j){

			
			compare[j][i]=AlignPWMs(allPWMs[i],allPWMs[j]);
			
			

		}
		

	}

}

double computeDr(Cluster &cluster,std::vector<std::vector<String<double> > > copyCompare){

	double Dr=0;
	
	for(unsigned i=0;i<length(cluster.content);++i){
		
		for(unsigned j=0;j<length(cluster.content);++j){
			
			
			if(cluster.content[i]<cluster.content[j]){
				

				Dr+= copyCompare[cluster.content[j]][cluster.content[i]][0];//copyCompare[i][j] --> i>j, daher falls i<j umdrehen
				
			}
			else if(cluster.content[j]<cluster.content[i]){
				
			
				Dr+= copyCompare[cluster.content[i]][cluster.content[j]][0];
				
			}
		}
	}
	


	return Dr/2;
}

double computeWk(int n,String<Cluster> &cluster){// n=allPWMsLength-n

	/***

		Wenn ein neues Cluster hinzukommt einfach addieren/n
		Wenn ein bestehendes verändert wird, den wert austauschen
		--> anhand von left/right machbar --> falls dort negative zahlen, dann verändert 
			falls nicht addieren 
			falls negativ --> -negativ -1 bestimmt cluster-nummer, dort Dr wert


	***/


	if(n>0){
		cluster[n].Wk=cluster[n-1].Wk;
		std::cout<<"Wk "<<cluster[n].Wk<<" "<<n<<" "<<cluster[n].left<<" "<<cluster[n].right<<" "<<cluster[n].Dr<<std::endl;
	}
	else
		cluster[n].Wk=0;

	if(cluster[n].left>=0 && cluster[n].right>=0)
		cluster[n].Wk +=cluster[n].Dr/length(cluster[n].content);//braucht nicht durch 2 dividieren, da in content nur die eine hälfte aufaddiert wurde
	else{
		if(cluster[n].left<0)
			cluster[n].Wk=cluster[n].Wk - cluster[-cluster[n].left-1].Dr/length(cluster[-cluster[n].left-1].content);//das Wk des vorherigen clusters abziehen --> 12 wird mit 3 zu 123 geklustert-->wäre 123.Dr/1 --> 12 wird von 12 und 3 abgezogen und 123 am ende der else addiert
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
			Für jede Spalte eines PWMs Normalverteilte Zufallszahlen
			So viele PWMs wie auch gemessen

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
		//std::cout<<"i "<<i<<std::endl;
		for(unsigned j=0;j<PWMLength;++j){
			ACGT="ACGT";
			uniDouble=0;
			for(unsigned k=0;k<3;++k){
				
				Pdf<Uniform<int> > uniformInt(0, length(ACGT)-1);
			
				Pdf<Uniform<double> > uniformDouble(0, 0.99-uniDouble);
				pos=pickRandomNumber(rng, uniformInt);
				ReferenzFreq[j][ACGT[pos]]=pickRandomNumber(rng, uniformDouble);
				uniDouble+=ReferenzFreq[j][ACGT[pos]];
				//std::cout<<"Pos: "<<ACGT[pos]<<" "<<uniDouble<<" "<<ReferenzFreq[j][ACGT[pos]]<<std::endl;
				erase(ACGT,pos);
			}
			ReferenzFreq[j][ACGT[0]]=1-uniDouble;
			//std::cout<<"Pos: "<<ACGT[0]<<" "<<uniDouble<<" "<<ReferenzFreq[j][ACGT[0]]<<std::endl;
			clear(ACGT);
			//std::cout<<std::endl;
	//		ReferenzFreq[j]["zufallACGT"]="zufallvon 0 bis 0.99";
	//		ReferenzFreq[j]["zufallACG"]="zufallvon 0 bis 1-[zufallACGT]";
	//		ReferenzFreq[j]["zufallAC"]="zufallvon 0 bis 1-[zufallACGT]-[zufallACG]";
	//		ReferenzFreq[j]["A"]="zufallvon 0 bis 1-[zufallACGT]-[zufallACG]-[zufallAC]";

		} 
		
		appendValue(Reference,ReferenzFreq,Exact());
		
		clear(ReferenzFreq);
	}



}

void compute_l_quer(String<Cluster> &cluster, 
					String<double> &l_quer_for_k,
					String<String<double> > &allWk,
					unsigned i,
					int B,
					unsigned allPWMsLength){

		unsigned k=0;
		
		for(unsigned j=0;j<allPWMsLength-1;++j){//addiert die log(Wk's) aller B-Referenz-Daten auf --> nur die Summe wird benötigt
			/*****
					j != anzahl der cluster -->allPWMsLength - (j+1) ist Anzahl der Cluster
					Bei 4 Objekten entspricht k=0 dem Fall, dass 2 zusammengeclustert wurden und der Rest allein ist
					k=allPWMsLength-2 entspricht einem großem Cluster
			*****/
			k=allPWMsLength-(j+1);
			appendValue(allWk[i],cluster[j].Wk);
			
			l_quer_for_k[0]+=log(cluster[j].Wk)/B;
			std::cout<<l_quer_for_k[0]<<std::endl;

			if(i==0){
				l_quer_for_k[k]=log(cluster[j].Wk)/B;
			}
			else{
				l_quer_for_k[k]+=log(cluster[j].Wk)/B;
				std::cout<<"else "<<k<<" "<<l_quer_for_k[k]<<" "<<log(cluster[j].Wk)<<std::endl<<std::endl;
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
		std::cout<<l_quer_for_k[k]<<" "<<log(observedCluster[j].Wk)<<" ";
		Gap[k]=l_quer_for_k[k]-log(observedCluster[j].Wk);
		std::cout<<Gap[k]<<std::endl;
	}

}

void compute_sdk_and_sk(String<double> &sk,
						String<double> &l_quer_for_k,
						String<String<double> > &allWk,
						unsigned allPWMsLength,
						int B){

	String<double> sdk;//Standard-Abweichungen
	resize(sdk,allPWMsLength);

	for(unsigned i=0;i<B;++i){

		for(unsigned k=0;k<allPWMsLength-1;++k){
			
			if(i==0)
				sdk[k]=(log(allWk[i][k])-l_quer_for_k[allPWMsLength-1])*(log(allWk[i][k])-l_quer_for_k[allPWMsLength-1]);
			else
				sdk[k]+=(log(allWk[i][k])-l_quer_for_k[allPWMsLength-1])*(log(allWk[i][k])-l_quer_for_k[allPWMsLength-1]);
		}
	}
	

	/*******************
		sk berechnen
	********************/

	double x=1+1/B;
	for(unsigned k=0;k<allPWMsLength-1;++k){

		sdk[k]=sqrt(sdk[k]/B);
		sk[k]=sdk[k]*sqrt(x);

	}
	clear(sdk);

}

void computeGapStat(String< std::map<unsigned int,std::map<Iupac,double> > > &allPWMs,
					String<Cluster> &observedCluster){

	String<Cluster> cluster;
	String<String<double> > allWk;
	String<double> Gap;
	String<double> l_quer_for_k;
	String<double> sk;
	int B=5;//B Referenz-Datensätze werden erstellt
	

	String< std::map<unsigned int,std::map<Iupac,double> > > Reference;
	unsigned allPWMsLength=length(allPWMs);
	unsigned PWMLength=length(allPWMs[0]);
	
	resize(allWk,B);
	
	resize(sk,allPWMsLength);
	resize(Gap,allPWMsLength);
	resize(l_quer_for_k,allPWMsLength);
	l_quer_for_k[0]=0;
	
	/************
		l_quer berechnen
	************/
	for(unsigned i=0;i<B;++i){

		computeReferenceData( Reference,allPWMsLength, PWMLength );
		
		PWMClustering(Reference,cluster);
		compute_l_quer(cluster,l_quer_for_k,allWk,i,B,allPWMsLength);
		clear(Reference);
		clear(cluster);
	}
	
	/******************
		Gap berechnen
	******************/
	computeGap(l_quer_for_k,observedCluster,Gap,allPWMsLength);


	/******************
		sdk && sk berechnen
	******************/

	compute_sdk_and_sk(sk,l_quer_for_k,allWk,allPWMsLength,B);

	/*******************
		kleinstes k für die Gap berechnen
	********************/
	for(unsigned k=0;k<allPWMsLength-2;++k){
		std::cout<<"gap "<<Gap[k]<<" "<<Gap[k+1]<<" "<<sk[k+1]<<std::endl;
		if(-Gap[k]>=-Gap[k+1]-sk[k+1]){
			std::cout<<"k "<<k<<std::endl;
			//break;
		}

	}


	
	clear(allWk);
	clear(Gap);
	



}

void PWMClustering(String< std::map<unsigned int,std::map<Iupac,double> > > &allPWMs,
				   String<Cluster> &cluster){
	
	unsigned allPWMsLength=length(allPWMs);	
	String<double> minDifference;
	String<int> traceback;
	std::vector<std::vector<String<double> > > compare(allPWMsLength, std::vector<String<double> >(allPWMsLength));
	std::vector<int> clusterId(allPWMsLength);

	resize(cluster,allPWMsLength);
	resize(minDifference,3);
	minDifference[0]=0;
	 
	
	unsigned j;
	for (j = 0; j < allPWMsLength; j++) clusterId[j] = j;//to assign which PWM is in which cluster
	/***
			Average Linkage:
	***/
	std::vector<unsigned> weights(allPWMsLength,1);
	
	


	//vorher per local alignment feststellen wo die beste Überlappung ist!
	computesDistantMatrix(allPWMs,compare,allPWMsLength);

	std::vector<std::vector<String<double> > > copyCompare(compare);
	String<int> temp;
	minDifferenceInMatrix(allPWMsLength,minDifference,compare);
	for(unsigned n=allPWMsLength;n>1;--n){//treshold noch bestimmen --> && minDifference[0]<2
		
		
		//BuildMeanOf2PWMs(seq,seq.allPWMs[int(minDifference[1])],seq.allPWMs[int(minDifference[2])]);//bildet aus 2PWMs die Mittelwerte und speichert sie in
		
		UpdateDistantMatrix(n,int(minDifference[1]),int(minDifference[2]),compare,weights,AverageLinkage());
	
		
		/****
				Erzeuge Vector mit allen Objekten des Clusters
				-> für computeDr
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
					Dr berechnet
		****/
		cluster[allPWMsLength-n].left=clusterId[int(minDifference[1])];
		cluster[allPWMsLength-n].right=clusterId[int(minDifference[2])];
		cluster[allPWMsLength-n].Dr=computeDr(cluster[allPWMsLength-n],copyCompare);
		
		/****
					Wk berechnet
		****/
		
		computeWk(allPWMsLength-n,cluster);
		

		appendValue(traceback,clusterId[int(minDifference[1])]);
		appendValue(traceback,clusterId[int(minDifference[2])]);
		/****
			Id's neu berechnen
		****/
		clusterId[int(minDifference[1])]=n-allPWMsLength-1;
		for(j=0;j+minDifference[2]<allPWMsLength-1;++j)
			clusterId[int(minDifference[2])+j]=clusterId[int(minDifference[2])+j+1];
		
		minDifferenceInMatrix(n-1,minDifference,compare);
	}
	
	clear(clusterId);
	clear(minDifference);
	clear(compare);
	clear(copyCompare);

	Iterator<String<int> >::Type tracebackIt;
	std::cout<<"Traceback: "<<std::endl;
	for(tracebackIt=begin(traceback);tracebackIt!=end(traceback);++tracebackIt){

		std::cout<<*tracebackIt<<" ";
		++tracebackIt;
		std::cout<<*tracebackIt<<std::endl;
	}

}


/*Prints the Mapping:
Kmer	Seq1	Seq2	...	Seqn	CumulatedCounter
-->Template

*/
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
	//std::cout<<pValueMap.size()<<std::endl;
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


//Test the Map-lengths match eachother and with the sequences
void DebugMap(  Seq &seq,
				Seq &back,
				std::map<String<Dna5>,std::vector<int> > &sequencesCounter,
				std::map<String<Dna5>,std::vector<int> > &backgroundCounter){

	typedef std::map<String<Dna5>,std::vector<int> > Dna5CounterMap;
	Dna5CounterMap::iterator MapIterator;
	MapIterator=sequencesCounter.begin();
	Dna5CounterMap::iterator MapIteratorB;
	MapIteratorB=backgroundCounter.begin();

	
	SEQAN_ASSERT_EQ(length(sequencesCounter),length(backgroundCounter));
	SEQAN_ASSERT_EQ(length((*MapIterator).second),(length(seq.ids)+1));//+1, because of the last field in vector 
	SEQAN_ASSERT_EQ(length((*MapIteratorB).second),(length(back.ids)+1));

	//std::cout<<length(sequencesCounter)<<std::endl;
	//std::cout<<length((*MapIterator).second)<<std::endl;
	//std::cout<<length(backgroundCounter)<<std::endl;
	//std::cout<<length((*MapIteratorB).second)<<std::endl;
	//std::cout<<length(seq.ids)<<std::endl;
	//std::cout<<length(back.ids)<<std::endl;
}

void DebugMultiMap( std::map<String<Dna5>,std::vector<int> > &sequencesCounter,
					std::multimap<double,String<Dna5> > &SortedPValue){
	SEQAN_ASSERT_EQ(length(sequencesCounter),SortedPValue.size());

}

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

double calcFET( unsigned int a,
				unsigned int b,
				unsigned int c,
				unsigned int d){

	return exp(logFac[a+b] + logFac[c+d] + logFac[a+c] + logFac[b+d] - 
		(logFac[a+b+c+d] + logFac[a]+ logFac[b] + logFac[c] +logFac[d]));
}


void modifyFET( unsigned int a,
				unsigned int b,
				unsigned int c,
				unsigned int d, 
				double &pValue){

	
	
	pValue= calcFET(a,b,c,d);
	//	//std::cout<<(*MapI).first<<" "<<pValue<<"   ";
	
		while(b!=0 && c!=0){//modify to be more extrem
			++a;
			--b;
			--c;
			++d;
			pValue += calcFET(a,b,c,d);
			
	//		//std::cout<<pValue<<"   ";
		}



}




/***************************

log((a+b)!(c+d)!(a+c)!(b+d)!/a!b!c!d!n!) = logFactorial(a+b) + logFactorial(c+d) +
										   logFactorial(a+c) + logFactorial(b+d) -
										   (logFactorial(a+b+c+d) + logFactorial(a)+
										   logFactorial(b) + logFactorial(c) +
										   logFactorial(d))

pValue = exp(log((a+b)!(c+d)!(a+c)!(b+d)!/a!b!c!d!n!))

a=In Sequenz gefunden	b=In Background gefunden
c=In Sequenz !gefunden	d=In Background !gefunden

a = sequenceCounter		b= backgroundCounter
a+c=SeqsNumber			b+d=backgroundNumber
--> c= SeqsNumber - cumulated(sequenceCounter)
	d= backgroundNumber - cumulated(backgroundCounter)

F√ºr den einseitigen Test zus√§tzlich:
	++a und --c
F√ºr den zweiseitigen Test:
	++a und --c
	--a und ++c



****************************/


void FisherExactTest(Seq &seq, 
					 Seq &back){

	
	

	double pValue=0;
	typedef std::map<String<Dna5>,unsigned int >::iterator MapIterator;
	MapIterator MapI=seq.seqCounter.begin();
	MapIterator MapIB=back.seqCounter.begin();
	//std::cout<<(*MapI).first<<" "<<(*MapI).second.back()<<std::endl;
	//std::cout<<(*MapIB).first<<" "<<(*MapIB).second.back()<<std::endl;
	
	for(;MapI!=seq.seqCounter.end();++MapI,++MapIB){
		

		modifyFET((*MapI).second,(*MapIB).second,(seq.SeqsNumber - (*MapI).second),(back.SeqsNumber - (*MapIB).second),pValue);
	
	
		//std::cout<<pValue<<std::endl;
		//SortedPValue[pValue]=(*MapI).first;
		seq.SortedPValue.insert(std::pair<double,String<Dna5> > (pValue, (*MapI).first));
		seq.SortedPValueReversed.insert(std::pair<String<Iupac>,double > ((*MapI).first,pValue));
	}


}
//--> templates verwenden
void FisherExactTest(std::map<String<Iupac>,unsigned int > &SequenceCounter,
					 std::map<String<Iupac>,unsigned int > &BackgroundCounter,
					 Seq &seq, 
					 Seq &back){

	
	

	double pValue=0;
	typedef std::map<String<Iupac>,unsigned int >::iterator MapIterator;
	MapIterator MapI=SequenceCounter.begin();
	MapIterator MapIB=BackgroundCounter.begin();
	//std::cout<<(*MapI).first<<" "<<(*MapI).second.back()<<std::endl;
	//std::cout<<(*MapIB).first<<" "<<(*MapIB).second.back()<<std::endl;
	
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
	
	
	

	//std::cout<<"begin Fisher "<<seq.generalizedKmer.size()<<" "<<back.generalizedKmer.size()<<std::endl;
	
	if(seq.generalizedKmer.size()==0)	
		return 2;

	typedef std::map<String<Iupac>,unsigned int >::iterator MapIterator;
	
	std::multimap<double,String<Iupac> >::iterator MapIterator2;
	MapIterator MapI = seq.generalizedKmer.begin();
	MapIterator MapIB= back.generalizedKmer.begin();
	double pValue=0;
	//std::cout<<generalizedKmerSequence.size();
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
	std::multimap<double,String<Dna5> >::iterator MapIterator;
	std::multimap<double,String<Iupac> >::iterator MapIteratorT;	
	std::multimap<double,String<Iupac> > generalizedSortedPValueTemp;
	std::map<String<Iupac>,unsigned int> generalizedKmerTemp;
	std::map<String<Iupac>,unsigned int> generalizedKmerBackgroundTemp;
	unsigned int i=0;
	unsigned int limit;
	if(seq.SortedPValue.size()>seq.seed)	limit=seq.seed;//seed meist = 100
	else if(seq.SortedPValue.size()==0) return;
	else	limit = seq.SortedPValue.size();
	for(MapIterator=seq.SortedPValue.begin();i<limit;++MapIterator,++i){//iterate over Top100
		GeneralizeKmer((*MapIterator).second,IMaps,seq,back);
	}
	/*
		- only do the next function call, if in the last at least one pValue<treshold	
		- call GeneralizeKmer in loop
	*/
	//PrintMap(seq.generalizedKmer);
	//seq.generalizedSortedPValue.insert(seq.SortedPValue.begin(),seq.SortedPValue.end());
	double topPValue = FisherExactTest(seq,back,generalizedSortedPValueTemp);// lowest pValue from the first generalization
	double topPValueOld =seq.SortedPValue.begin()->first;//lowest pValue before generalization
	
	while(topPValue<0.05 && topPValue<topPValueOld){//only start a new round, if the top PValue is an improvement of the old one
		
		
		/*
			while wird das erste mal mit generalizedKmer aufgerufen und dem tempor√§ren mapping der pValues
			das tempor√§re mapping wird in das richtige mapping gemerged und gecleant, damit geschaut werden kann, ob bei den neuen pValues ein wert
			√ºber dem treshold ist --> falls nicht bricht die while ab
			falls doch wird generalizedKmer kopiert und gecleant aus dem gleichen grund, damit, nur die neuen generalisierten Kmere untersucht werden
			--> das Temp hier, um √ºber alle alten zu gehen, um diese weiter zu generalisieren

		*/

		seq.generalizedSortedPValue.insert(generalizedSortedPValueTemp.begin(),generalizedSortedPValueTemp.end());

	
		
	
		generalizedKmerTemp.clear();
		generalizedKmerBackgroundTemp.clear();
		generalizedKmerTemp=seq.generalizedKmer;
		generalizedKmerBackgroundTemp=back.generalizedKmer;
		back.generalizedKmer.clear();
		seq.generalizedKmer.clear();
	//	generalizedMapIterator= generalizedKmerTemp.begin();
		if(generalizedSortedPValueTemp.size()>seq.seed)	limit=seq.seed;
		else if(generalizedSortedPValueTemp.size()==0) return;
		else	limit = generalizedSortedPValueTemp.size();
	
		i=0;//only Top100
		for(MapIteratorT=generalizedSortedPValueTemp.begin();i<limit;++MapIteratorT,++i){//iterate over Top100
		
			//Temp ums zu finden, aber das normale auch √ºbergeben, zum neu bef√ºllen
			
			GeneralizeKmer((*MapIteratorT).second,generalizedKmerTemp,generalizedKmerBackgroundTemp,IMaps,seq,back);
			
			                                                                        
	
		}
		generalizedSortedPValueTemp.clear();
		//std::cout<<"nach for"<<std::endl;
	topPValueOld =topPValue;
	topPValue = FisherExactTest(seq,back,generalizedSortedPValueTemp);
	//std::cout<<"nach Fisher "<<topPValue<<" "<<topPValueOld<<std::endl;
		
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
			//std::cout<<temp<<" ";
			tempChar =*tempIt;//stores the current char because of temp2
			*tempIt=*replaceIt;
			temp2=temp;//temp2 ist nun das mit dem char für die neue wildcard
			*tempIt=tempChar;//temp wieder das alte, wird aber im nächsten schritt mit einer neuen wildcard ergänzt
			*tempIt =IMaps.IupacMap[IMaps.IupacMapReversed[*tempIt] + IMaps.IupacMapReversed[*replaceIt]];//compute Iupac-letter--> A + G = R and replace the current location in temp
			
			//std::cout<<Kmer<<" "<<temp<<" "<<temp2<<std::endl;
			if(seq.generalizedKmer.find(temp)!=seq.generalizedKmer.end()) continue;// if Kmer is in the Map -->nothing to do
			//estimateCounter mit Kmer und temp2 aufrufen --> Kmer=AAA  temp2=TAA temp=WAA
			estimateCounter(seq,Kmer,temp2,counter);
			
			seq.generalizedKmer[temp]=counter;//temp ist das neue motiv
			estimateCounter(back,Kmer,temp2,counter);
			back.generalizedKmer[temp]=counter;
		}
	}

	//PrintMap(seq.generalizedKmer);
	

}

/*	- the same as above except that each String has already a wildcard
*/
void GeneralizeKmer(String<Iupac> Kmer,
					std::map<String<Iupac>,unsigned int> &generalizedKmerTemp,
					std::map<String<Iupac>,unsigned int> &generalizedKmerBackgroundTemp,
					IupacMaps &IMaps,
					Seq &seq, 
					Seq &back){
	String<Iupac> temp;//temporary String --> generalizedKmer[temp]
	String<Iupac> temp2;//Kmer with replaced position--> relevant for estimateCounter
	Iterator<String<Iupac> >::Type tempIt;//Iterator over temp --> same length as Kmer
	String<Iupac> replace;
	Iterator<String<Iupac> >::Type replaceIt;
	unsigned int counter =0;
	char tempChar;
	
	temp = Kmer;
	tempIt = begin(temp);
	//std::cout<<temp<<" ";
	for(;tempIt!=end(temp);++tempIt){//loop over each position in kmer
		//if(*tempIt == 'A' || *tempIt == 'C' ||*tempIt == 'G' ||*tempIt == 'T') continue;//only replace the position with a wildcard
		
		
		if(*tempIt =='N')continue;//gibt nichts mehr zu ersetzen
		replace=IMaps.IupacMapReplace[*tempIt];
		replaceIt = begin(replace);	
		for(;replaceIt!=end(replace);++replaceIt){// loop over the replacement-chars, W=AT --> replace=CG
			temp = Kmer;// reset temp
			//std::cout<<*replaceIt<<std::endl;
			tempChar =*tempIt;//stores the current char because of temp2
			*tempIt=*replaceIt;//replace the current char for temp2
			temp2=temp;
			*tempIt=tempChar;

			*tempIt =IMaps.IupacMap[IMaps.IupacMapReversed[*tempIt] + IMaps.IupacMapReversed[*replaceIt]];
			
			if(seq.SortedPValueReversed[Kmer] >= 0.05 || seq.SortedPValueReversed[temp2] >= 0.05) continue;//only if Kmer and temp2 are significant estimate the counter
			if(generalizedKmerTemp.find(temp)!=generalizedKmerTemp.end()) continue;
			//std::cout<<Kmer<<" "<<temp<<" "<<temp2<<std::endl;
			estimateCounter(seq,generalizedKmerTemp,Kmer,temp2,counter);
			seq.generalizedKmer[temp]=counter;
			estimateCounter(back,generalizedKmerBackgroundTemp,Kmer,temp2,counter);
			back.generalizedKmer[temp]=counter;
			//std::cout<<temp<<" "<<counter<<std::endl;
		}

		}
		/**
		- gehe √ºber Kmer --> if ACGT continue
		- else rufe IupacMapReplace auf
		- -->for(;replaceIt!=end(replace);++replaceIt) --> loop √ºber den String aus IupacMapReplace
		- temp ist das neue Kmer und temp2 die neue wildcard--> estimateCounter aufrufen --> f√ºr fore- und background

		- funktion wird in for-schleife aufgerufen --> geht √ºber alle
		- clear generaliedKmer vor der for-schleife
		- --> wird neu bef√ºllt und kann mit FisherExact aufgerufen werden
		- SortedPValue wird berechnet (in temp) und √ºberpr√ºft ob noch ein wert drunter ist-->falls ja insertiere temp in den rest
		- falls nein nehme die top100 und rufe exactSearch auf




	**/




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
	//std::cout<<temp<<"  "<<temp2<<std::endl;
	unsigned int RE1=seq.seqCounter.find(temp)->second;//das alte motiv aus dem letzten schritt
	
	if(seq.seqCounter.find(temp2)!=seq.seqCounter.end()){//falls temp2 ein altes motiv ist, hat es einen counter
		counter= RE1+ seq.seqCounter.find(temp2)->second - (RE1*seq.seqCounter.find(temp2)->second)/seq.SeqsNumber;
			//std::cout<<temp<<" "<<temp2<<" "<<RE1<<"  "<<SequenceCounter.find(temp2)->second.back()<<" "<<counter<<std::endl;
		if(counter>seq.SeqsNumber){ std::cout<<"if "<<counter<<" "<<temp<<" "<<RE1<<" "<<temp2<<" "<<seq.seqCounter.find(temp2)->second<<" SeqsNumer "<<seq.SeqsNumber<<std::endl;
			system("PAUSE");
		}
	}
	else{
		counter=RE1;//RE2=0, da noch nicht vorhanden
		//std::cout<<"else "<<counter<<std::endl;
		if(counter>seq.SeqsNumber){ std::cout<<"else "<<counter<<" "<<temp<<" "<<RE1<<" "<<temp2<<" SeqsNumer "<<seq.SeqsNumber<<std::endl;
			system("PAUSE");
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

	//String<Iupac>="AAWR";
	
	unsigned int RE1=(*generalizedKmer.find(temp)).second;//the new seed RE is a Kmer with wildcard
	//temp2 may be in generalizedKmer or in SequenceCounter
	if(seq.seqCounter.find(temp2)!=seq.seqCounter.end()){// if temp2 is in SequenceCounter do the same as above --> has no wildcard
		counter= RE1+ seq.seqCounter.find(temp2)->second- (RE1*seq.seqCounter.find(temp2)->second)/seq.SeqsNumber;
		//std::cout<<temp<<" "<<temp2<<" "<<RE1<<"  "<<SequenceCounter.find(temp2)->second.back()<<" "<<counter<<std::endl;
		if(counter>seq.SeqsNumber){ std::cout<<"if "<<counter<<" "<<temp<<" "<<RE1<<" "<<temp2<<" "<<seq.seqCounter.find(temp2)->second<<" SeqsNumer "<<seq.SeqsNumber<<std::endl;
			system("PAUSE");
		}
	}
		
	else if(generalizedKmer.find(temp2)!=generalizedKmer.end()){//if temp2 has a wildcard and is found in generalizedKmer
		counter= RE1+ generalizedKmer.find(temp2)->second - (RE1*generalizedKmer.find(temp2)->second)/seq.SeqsNumber;
		//std::cout<<temp<<" "<<temp2<<" "<<RE1<<"  "<<generalizedKmer.find(temp2)->second<<" "<<counter<<std::endl;
		if(counter>seq.SeqsNumber){ std::cout<<"elif "<<counter<<" "<<temp<<" "<<RE1<<" "<<temp2<<" "<<generalizedKmer.find(temp2)->second<<" SeqsNumer "<<seq.SeqsNumber<<std::endl;
			system("PAUSE");
		}
		
	}
	else{//if temp2 is not found
		counter= RE1;//RE2=0
		if(counter>seq.SeqsNumber){ std::cout<<"else "<<counter<<" "<<temp<<" "<<RE1<<" "<<temp2<<" SeqsNumer "<<seq.SeqsNumber<<std::endl;
				system("PAUSE");
		}
	}
}


void exactGeneralizeCount(  std::map<String<Iupac>,unsigned int > &seqCounter,
							std::map<String<Iupac>,unsigned int > &backCounter,
							Finder<Index<StringSet<String<Dna5> > > > &finder,
							Finder<Index<StringSet<String<Dna5> > > > &finderB,
							Seq &seq,
							Seq &back,
							IupacMaps &IMap){

	std::multimap<double,String<Iupac> >::iterator generalizedSortedPValueIt;
	//std::map<String<Iupac>,double > generalizedSortedPValueReversed;

	generalizedSortedPValueIt = seq.generalizedSortedPValue.begin();

	for(unsigned int i=0;i<seq.seed && generalizedSortedPValueIt!=seq.generalizedSortedPValue.end() ;++i,++generalizedSortedPValueIt){
		//std::cout<<length(seq.generalizedSortedPValue)<<" "<<seq.generalizedSortedPValue.size()<<" ";
		//std::cout<<(*generalizedSortedPValueIt).second<<" ";
		if(seqCounter.find((*generalizedSortedPValueIt).second)!=seqCounter.end()) continue;
		
		CountKmer(seqCounter,finder,(*generalizedSortedPValueIt).second,seq,IMap);
		CountKmer(backCounter,finderB,(*generalizedSortedPValueIt).second,back,IMap);
		
	}
	
	seq.generalizedSortedPValue.clear();
	
	//PrintMap(seqCounter,seq.SeqsNumber);
	
	FisherExactTest(seqCounter,backCounter,seq,back);//computes the pValue of each Motif due to the counter
	std::cout<<std::endl;
	//PrintMap(seq.generalizedSortedPValue);
	seqCounter.clear();
	backCounter.clear();



}



String<double> AlignPWMs(std::map<unsigned int,std::map<Iupac,double> > &freqMatrix1,std::map<unsigned int,std::map<Iupac,double> > &freqMatrix2){

  int freqL1 = length(freqMatrix1);                     
  int freqL2 = length(freqMatrix2);

  std::vector<std::vector<double> > M(freqL1+1,std::vector<double>(freqL2+1));     
  for(unsigned i=0;i<=freqL1;++i){
    for(unsigned j=0;j<=freqL2;++j){
      M[i][j]=0;
    }
  } 
 
  
	//traceback wird nicht benötigt, da zwischen den anfangs- und end-gaps keine erlaubt sind
	String<double> Mmax;
	appendValue(Mmax,100);
	appendValue(Mmax,0);
	appendValue(Mmax,0);

	
	for(unsigned i=1;i<=freqL1;++i){
		for(unsigned j=1;j<=freqL2;++j){

			M[i][j]=M[i-1][j-1]+ComparePWM(freqMatrix1[i-1],freqMatrix2[j-1],Entropy());//je größer, desto unterschiedlicher --> Problem: je weniger abgezogen wird, desto größer
			
			
				
			if((i==freqL1 && j>freqL2*0.5  && M[i][j]/j<Mmax[0])){//auf die Länge der Überlappung normalisieren --> /j --> über die Hälfte des Kmers soll überlappen
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
