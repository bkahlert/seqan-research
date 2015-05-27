#ifndef SANDBOX_MEYERCLP_APPS_DREME_H_
#define SANDBOX_MEYERCLP_APPS_DREME_H_

#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <math.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/stream.h>
#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/misc/misc_interval_tree.h>

using namespace seqan;


struct Seq
{
	StringSet<CharString> ids;
	unsigned c; // Anzahl der Motive, provisorisch
	StringSet<String<Dna5> > seqs;//
	Index< StringSet<String<Dna5> > > SArray;
	
	unsigned int SeqsNumber;
	
	std::map<String<Dna5>,unsigned int > seqCounter;//maps the Sequence-Kmere to a Counter for the Sequence
	std::map<String<Iupac>,unsigned int> generalizedKmer;

	std::map<String<Iupac>,double > SortedPValueReversed;
	
	std::multimap<double,String<Dna5> > SortedPValue;
	std::multimap<double,String<Iupac> > generalizedSortedPValue;

	/**
		In BuildFrequencyMatrix werden die gefundenen Intervalle des Top-Motivs an Intervalls gehängt
		Am Ende der while wird dann ein IntervallTree für jede Sequenz erzeugt, sodass im nächsten Schritt direkt geschaut werden kann ob ein Motiv mit einem schon 
		gefundenen überlappt
	**/
	typedef IntervalAndCargo<unsigned, unsigned> TInterval;
	String<String<TInterval> > intervals;//String<TInterval> enthält alle Intervalle einer Sequenz. String<String<..>> enthält alle Sequenzen
	
	String<IntervalTree<unsigned> > intervalTrees; // Tree für schnellere Suche erstellen --> String von Trees, da mehrere Sequenzen
	String<unsigned> results;
};

struct IupacMaps
{
	std::map<unsigned int,char> IupacMap;
	std::map<char,unsigned int> IupacMapReversed;
	std::map<char,String<Iupac> > IupacMapReplace; //stores the replacement-chars
	std::map<char,String<Dna5> > IupacMapReplaceReversed;
	std::map<String<Dna5>, char > IupacMapInversed;
};


void readFastA(struct Seq &seq, CharString fname);
template <typename TStream>
void PrintFastA(TStream & stream, Seq &seq);
void initST(Seq &seq);
void PrintST(Seq &seq);
void initExactKmer(Seq &seq,Seq &back,unsigned int kmer_len,unsigned int kmer_len_end);
void CountKmer(Seq &seq, Finder<Index<StringSet<String<Dna5> > > > &finder, String<Dna5> &Kmer);
void CountKmer(std::map<String<Iupac>,unsigned int > &Dna5CounterMap, Finder<Index<StringSet<String<Iupac> > > > &finder, String<Iupac> &Kmer,Seq &seq,IupacMaps &IMap);
void PrintMap(std::map<String<Dna5>,unsigned int > &Dna5CounterMap,unsigned int SeqsNumber);
void PrintMap(std::map<String<Iupac>,unsigned int > &Dna5CounterMap,unsigned int SeqsNumber);
void PrintMap(std::multimap<double,String<Dna5> > &pValueMap);
void PrintMap(std::map<String<Iupac>,unsigned int> &generalizedKmer);
void PrintMap(std::map<unsigned int,std::map<Iupac,double> > freqMatrix,unsigned int Kmerlength,bool foreground);
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
void FindTopKmer(Seq &seq,String<Iupac> &temp,Finder<Index<StringSet<String<Dna5> > > > &finder,unsigned int &counter,std::vector<int> &CounterV,String<unsigned int> &replaceString,std::map<unsigned int,std::map<Iupac,double> > &freqMatrix);
void loopOverTopKmer(	Seq &seq,String<Iupac> &temp,String<Iupac> &Kmer,Iterator<String<Iupac> >::Type &tempIt,Finder<Index<StringSet<String<Dna5> > > > &finder,unsigned int &counter,std::vector<int> &CounterV,IupacMaps &IMap,std::map<unsigned int,std::map<Iupac,double> > &freqMatrix,String<unsigned int> &replaceString);
void BuildFrequencyMatrix(std::map<unsigned int,std::map<Iupac,double> > &freqMatrix,std::map<String<Iupac>,unsigned int > &seqCounter, Finder<Index<StringSet<String<Dna5> > > > &finder, String<Iupac> &Kmer,Seq &seq, IupacMaps &IMaps,String<unsigned int> &replaceString);
void replaceKmer(Seq &seq,unsigned int stringNumber, unsigned int begin, unsigned int end);



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
				if(CounterV[beginPosition(finder).i1] == 0 || seq.results>0){//count number of sequences containing the motif, not the occurrences to avoid problems with self-overlapping
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
						if(CounterV[beginPosition(finder).i1] == 0 || seq.results>0){//count number of sequences containing the motif, not the occurrences to avoid problems with self-overlapping
							//ansonsten m√ºsste das array noch einmal durch gegangen werden und an jeder stellt !=0 ++
							++CounterV[beginPosition(finder).i1];
							++CounterV[seq.SeqsNumber];//last Position in CounterV is cumulated sum
							++counter;
							if(temp=="CGTW")
								std::cout<<temp<<" ";
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
			//std::cout<<temp<<" "<<CounterV[SeqsNumber]<<std::endl;
			CounterV.clear();
			
			/*
				Der clear hier bewirkt, dass bei ASG entweder ACG oder AGG vorkommen darf,
				falls beide vorkommen(in einer Sequenz) wird der counter trotzdem nur um 1 erh√∂ht
				-->Counter f√ºr AGG ist falsch, ASG stimmt jedoch wieder. Counter[AGG] --> irrelevant
				-->AGG wird jedes mal wenns ben√∂tigt wird neu berechnet-->optimieren
			*/
		

}



void replaceKmer(   Seq &seq,
					String<unsigned int> &replaceString){

	unsigned int beg =0;
	unsigned int endP =0;
	unsigned int Snumber =0;

	//std::cout<<stringNumber<<" "<<begin<<" "<<end<<std::endl;
	Iterator<String<unsigned int> >::Type StringIt;
	StringIt = begin(replaceString);
	for(;StringIt!=end(replaceString);++StringIt){
		Snumber = *StringIt;
		++StringIt;
		beg		= *StringIt;
		++StringIt;
		endP		= *StringIt;
		//std::cout<<Snumber<<" "<<beg<<" "<<endP<<std::endl;
		for(;beg<endP;++beg)
		{
			seq.seqs[Snumber][beg]='N';
		}
			//replace(seq.seqs[*StringIt],*(++StringIt),*(++StringIt),'N');
	}


}


///////////////////////////////////
void FindTopKmer(Seq &seq,
				String<Iupac> &temp,
				Finder<Index<StringSet<String<Dna5> > > > &finder,
				unsigned int &counter,
				std::vector<int> &CounterV,
				String<unsigned int> &replaceString,
				std::map<unsigned int,std::map<Iupac,double> > &freqMatrix){
	typedef IntervalAndCargo<unsigned, unsigned> TInterval;				
	clear(finder);
	//std::cout<<temp<<" vor-while"<<std::endl;
	while(find(finder,temp)){//search the current Kmer in all sequences
		//std::cout<<" "<<temp<<" "<<Kmer<<" "<<*replaceIt<<std::endl;
		//std::cout<<'[' <<beginPosition(finder)<<','<<endPosition(finder)<<")\t"<<infix(finder)<<std::endl;//Debug
		//replaceKmer(seq,beginPosition(finder).i1, beginPosition(finder).i2, endPosition(finder).i2);
		//std::cout<<temp<<" while"<<std::endl;
		//std::cout<<beginPosition(finder).i1<<" "<<beginPosition(finder).i2<<" "<<endPosition(finder).i1<<" "<<endPosition(finder).i2<<std::endl;
		//appendValue(replaceString,beginPosition(finder).i1);//i1 = Sequenz nummer, i2= stelle in der Sequenz
		//appendValue(replaceString,beginPosition(finder).i2);
		//appendValue(replaceString,endPosition(finder).i2);
		
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
		freqMatrix[k][temp[k]]+=CounterV[seq.SeqsNumber]; //GCAGCA --> counter der einzelnen wird um die gleiche anzahl hochgezählt, GCAGTA --> usw. Nicht-Wildcards haben W'keit 1
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
						IupacMaps &IMaps,
						std::map<unsigned int,std::map<Iupac,double> > &freqMatrix,
						String<unsigned int> &replaceString){

						
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
				loopOverTopKmer(seq,temp,temp,++tempIt,finder,counter,CounterV,IMaps,freqMatrix,replaceString);//only replace the position with a wildcard
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
						loopOverTopKmer(seq,temp,temp,++tempIttemp,finder,counter,CounterV,IMaps,freqMatrix,replaceString);
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
						FindTopKmer(seq,temp,finder,counter,CounterV,replaceString,freqMatrix);
						
					}
					*tempIt=resetTemp;
					//if ende von replaceIt, dann begin(replace) in temp speichern und mitübergeben--> referenz?
				}


}


void BuildFrequencyMatrix(  std::map<unsigned int,std::map<Iupac,double> > &freqMatrix,
							std::map<String<Iupac>,unsigned int > &seqCounter,
							Finder<Index<StringSet<String<Dna5> > > > &finder,
							String<Iupac> &Kmer,
							Seq &seq, 
							IupacMaps &IMaps,
							String<unsigned int> &replaceString){
			//std::cout<<Kmer<<std::endl;
			//freqMatrix -->unsigned int = position in Kmer, position 1 in map = prob. for A, pos. 2 = prob. for C...
			String<Iupac> temp;	
			
			Iterator<String<Iupac> >::Type tempIt;
			
			temp=Kmer;
			tempIt = begin(temp);
			unsigned int counter=0;
			
			std::vector<int> CounterV(seq.SeqsNumber+1,0);
		
			
			loopOverTopKmer(seq,temp,Kmer,tempIt,finder,counter,CounterV,IMaps,freqMatrix,replaceString);
			CounterV.clear();
			//loopOver funktionier, aber jetzt wird der counter nicht mehr richtig berechnet --> fixen  + andere loop anpassen
			//Durch die wildcards mehrere Vorkommen pro Sequence möglich:
			//seqCounter[temp]=CounterV[seq.SeqsNumber];
					
			//std::cout<<temp<<" "<<*replaceIt<<" ";
			
					
					
			
				
			if(counter>0){//normalisieren des Counters
				
				
				for( unsigned int k =0;k< length(temp);++k){
					
					freqMatrix[k]['A']=freqMatrix[k]['A']/counter; 
					freqMatrix[k]['C']=freqMatrix[k]['C']/counter; 
					freqMatrix[k]['G']=freqMatrix[k]['G']/counter; 
					freqMatrix[k]['T']=freqMatrix[k]['T']/counter; 
				}
	

				}
			else
				freqMatrix.clear();
				
				
	
			}
 			
}
//////////////////////////////////////

//FreqMatrix output
void PrintMap(  std::map<unsigned int,std::map<Iupac,double> > freqMatrix,
				unsigned int Kmerlength,
				bool foreground){
	std::map<Iupac,double> freq;
	std::cout<<std::endl;
	if(foreground)
		std::cout<<"foreground: "<<std::endl;
	else
		std::cout<<"background: "<<std::endl;
	for(unsigned int j=0;j<Kmerlength;++j){
				freq=freqMatrix[j];
				std::cout<<"Position: "<<j<<" A: "<<freq['A']<<std::endl;
				std::cout<<"Position: "<<j<<" C: "<<freq['C']<<std::endl;
				std::cout<<"Position: "<<j<<" G: "<<freq['G']<<std::endl;
				std::cout<<"Position: "<<j<<" T: "<<freq['T']<<std::endl;
				std::cout<<std::endl;
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
	
	
	

	//std::cout<<"begin"<<std::endl;
	if(seq.generalizedKmer.size()==0)	
		return 2;

	typedef std::map<String<Iupac>,unsigned int >::iterator MapIterator;
	
	std::multimap<double,String<Iupac> >::iterator MapIterator2;
	MapIterator MapI = seq.generalizedKmer.begin();
	MapIterator MapIB= back.generalizedKmer.begin();
	double pValue=0;
	//std::cout<<generalizedKmerSequence.size();
	for(;MapI!=seq.generalizedKmer.end();++MapI,++MapIB){
		

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
	if(seq.SortedPValue.size()>100)	limit=100;
	else	limit = seq.SortedPValue.size();
	for(MapIterator=seq.SortedPValue.begin();i<limit;++MapIterator,++i){//iterate over Top100
		GeneralizeKmer((*MapIterator).second,IMaps,seq,back);
	}
	/*
		- only do the next function call, if in the last at least one pValue<treshold	
		- call GeneralizeKmer in loop
	*/
	PrintMap(seq.generalizedKmer);
	//seq.generalizedSortedPValue.insert(seq.SortedPValue.begin(),seq.SortedPValue.end());
	double topPValue = FisherExactTest(seq,back,generalizedSortedPValueTemp);// lowest pValue from the first generalization
	double topPValueOld =seq.SortedPValue.begin()->first;//lowest pValue before generalization

	do{//only start a new round, if the top PValue is an improvement of the old one
		
		//std::cout<<"test";
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
		if(generalizedSortedPValueTemp.size()>100)	limit=100;
		else	limit = generalizedSortedPValueTemp.size();
	
		i=0;//only Top100
		for(MapIteratorT=generalizedSortedPValueTemp.begin();i<limit;++MapIteratorT,++i){//iterate over Top100
		
			//Temp ums zu finden, aber das normale auch √ºbergeben, zum neu bef√ºllen
			
			GeneralizeKmer((*MapIteratorT).second,generalizedKmerTemp,generalizedKmerBackgroundTemp,IMaps,seq,back);
			
			                                                                        
	
		}
		generalizedSortedPValueTemp.clear();

	topPValueOld =topPValue;
	topPValue = FisherExactTest(seq,back,generalizedSortedPValueTemp);
	
		
	}while(topPValue<0.05 && topPValue<topPValueOld);
	
	
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
			temp2=temp;
			*tempIt=tempChar;
			*tempIt =IMaps.IupacMap[IMaps.IupacMapReversed[*tempIt] + IMaps.IupacMapReversed[*replaceIt]];//compute Iupac-letter--> A + G = R and replace the current location in temp
			//std::cout<<Kmer<<" "<<temp<<" "<<temp2<<std::endl;
			if(seq.generalizedKmer.find(temp)!=seq.generalizedKmer.end()) continue;// if Kmer is in the Map -->nothing to do
			//estimateCounter mit Kmer und temp2 aufrufen
			estimateCounter(seq,Kmer,temp2,counter);
			//std::cout<<temp<<" ";
			seq.generalizedKmer[temp]=counter;
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
		
		
		
		replace=IMaps.IupacMapReplace[*tempIt];
		replaceIt = begin(replace);	
		for(;replaceIt!=end(replace);++replaceIt){// loop over the replacement-chars
			temp = Kmer;// reset temp
			//std::cout<<*replaceIt<<std::endl;
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
					 String<Iupac> temp,String<Iupac> temp2,
					 unsigned int &counter){
	//std::cout<<temp<<"  "<<temp2<<std::endl;
	unsigned int RE1=seq.seqCounter.find(temp)->second;
	unsigned int RE2=0;//counter for the second regular expression(Kmer)-->may be new
	if(seq.seqCounter.find(temp2)!=seq.seqCounter.end()){
		counter= RE1+ seq.seqCounter.find(temp2)->second - (RE1*seq.seqCounter.find(temp2)->second)/seq.SeqsNumber;
			//std::cout<<temp<<" "<<temp2<<" "<<RE1<<"  "<<SequenceCounter.find(temp2)->second.back()<<" "<<counter<<std::endl;
	}
	else{
		counter=RE1 + RE2;
		//std::cout<<"else "<<counter<<std::endl;
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
	unsigned int RE2=0;//temp2 may be in generalizedKmer or in SequenceCounter
	if(seq.seqCounter.find(temp2)!=seq.seqCounter.end()){// if temp2 is in SequenceCounter do the same as above --> has no wildcard
		counter= RE1+ seq.seqCounter.find(temp2)->second- (RE1*seq.seqCounter.find(temp2)->second)/seq.SeqsNumber;
		//std::cout<<temp<<" "<<temp2<<" "<<RE1<<"  "<<SequenceCounter.find(temp2)->second.back()<<" "<<counter<<std::endl;
		
	}
		
	else if(generalizedKmer.find(temp2)!=generalizedKmer.end()){//if temp2 has a wildcard and is found in generalizedKmer
		counter= RE1+ generalizedKmer.find(temp2)->second - (RE1*generalizedKmer.find(temp2)->second)/seq.SeqsNumber;
		//std::cout<<temp<<" "<<temp2<<" "<<RE1<<"  "<<generalizedKmer.find(temp2)->second<<" "<<counter<<std::endl;
	}
	else{//if temp2 is not found
		counter= RE1 + RE2;
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

	for(unsigned int i=0;i<100 && generalizedSortedPValueIt!=seq.generalizedSortedPValue.end() ;++i,++generalizedSortedPValueIt){
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

#endif  // #ifndef SANDBOX_MEYERCLP_APPS_DREME_H_