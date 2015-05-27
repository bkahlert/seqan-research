
#include "dreme.h"

using namespace seqan;



int main(int argc, char const ** argv){
	
	Seq sequences;
    Seq background;
	IupacMaps IMaps;
	unsigned int kmer_len=4;
	unsigned int kmer_len_end=6;
	sequences.seed=100;
	MapIupac(IMaps);//IupacMap for generalization
	
	if(argc !=3){
		std::cerr<<"ERROR: Invalid argument count."<<std:: endl
				 <<"Usage:" <<argv[0]<<"File"<<std::endl;
		return 1;
	}


	
	readFastA(sequences,argv[1]);
	//readFastA(background,argv[2]);
	//PrintFastA(sequences);//Debug
	background.seqs=sequences.seqs;
	sequences.SeqsNumber=length(sequences.seqs);//number of sequences
	background.SeqsNumber=length(background.seqs);//number of sequences
	
	for(unsigned int i=0;i<sequences.SeqsNumber-1;++i){
	
		std::random_shuffle(begin(background.seqs[i]),end(background.seqs[i]));
		std::random_shuffle(begin(background.seqs[i]),end(background.seqs[i]));
		std::random_shuffle(begin(background.seqs[i]),end(background.seqs[i]));
		std::random_shuffle(begin(background.seqs[i]),end(background.seqs[i]));
		//std::random_shuffle(begin(background.seqs[i]),end(background.seqs[i]));
	}
	
	/*std::ofstream outfile;
	outfile.open("test2.fasta");
	write2(outfile,sequences.ids,background.seqs,Fasta());

	outfile.close();
	*/
	std::ofstream output;
	


	
	//PrintFastA(background);//Debug
	logFactorial(sequences.SeqsNumber+background.SeqsNumber);//save all relevant factorial numbers
	resize(sequences.intervals, sequences.SeqsNumber);
	resize(background.intervals, background.SeqsNumber);
	resize(sequences.intervalTrees, sequences.SeqsNumber);
	resize(background.intervalTrees, background.SeqsNumber);

	sequences.c=1;
	background.c=1;
	initST(sequences);
	initST(background);


	/****
	Hintergrund-Wahrscheinlichkeit
	****/
	priorFreq(sequences);

	/****
	PSSM
	****/
	

	do{
		std::cout<<"new "<<std::endl;
	
	
		sequences.generalizedSortedPValue.clear();
		sequences.SortedPValue.clear();
		/*****
			- Bruteforce --> gets alls possible motifs form length 3 to 8 and counts them.
			- Only one occurence per sequence allowed
		*****/
		initExactKmer(sequences,background,kmer_len,kmer_len_end);

		std::cout<<"Exact done"<<std::endl;

		/*****
			- Computes the pValue of each motif due to the counter and saves it in SortedPValue
		*****/
	
		FisherExactTest(sequences,background);

		std::cout<<"Fisher done"<<std::endl;

		/*****
			- initiate by repeatedly picking the top motifs from SortedPValue
			- calls GeneralizeKmer and adds one wildcard to the Kmer and estimates the counter
			- after one generalization-round(top 100 motifs) the FisherExactTest is called, to compute the new SortedPValue-Map
			- if there is a pValue<treshold start again by picking the new top motifs
		*****/
		//PrintMap(sequences.SortedPValue);
		InitGeneralization(IMaps,sequences,background);
		std::cout<<"Generalize done"<<std::endl;
		//PrintMap(sequences.generalizedSortedPValue);
		//sequences.seqCounter.clear();
		//background.seqCounter.clear();
		sequences.SortedPValueReversed.clear();

		/*****
			- if there is not a single pValue<treshold exit the programm
		*****/
		std::map<unsigned int,std::map<Iupac,double> > freqMatrix;
		std::map<unsigned int,std::map<Iupac,double> > freqMatrixB;
		Finder<Index<StringSet<String<Dna5> > > > finder(sequences.SArray);
		Finder<Index<StringSet<String<Dna5> > > > finderB(background.SArray);//finder ins struct
		sequences.pValue=0;
		if(sequences.generalizedSortedPValue.size()>0){
			
	
		
			//PrintMap(sequences.generalizedSortedPValue);

			std::map<String<Iupac>,unsigned int > seqCounter;
			std::map<String<Iupac>,unsigned int > backCounter;
			
	
			/*****
				- gets the top 100(generalized) motifs and computes the exact counter and pValue
			*****/
	
			exactGeneralizeCount(seqCounter,backCounter, finder, finderB,sequences,background, IMaps);
			std::cout<<"exactGeneralize done"<<std::endl;

			//PrintFastA(sequences);

			//PrintMap(sequences.generalizedSortedPValue);

			//std::cout<<sequences.generalizedSortedPValue.begin()->second;
			output.open("output.fasta",std::ios::out|std::ios::app);
			write(output,(*sequences.generalizedSortedPValue.begin()).second);
			write(output,"   ");
		
			output<<(*sequences.generalizedSortedPValue.begin()).first;
			write(output,"   ");
			output.close();
			BuildFrequencyMatrix(freqMatrix, finder, (*sequences.generalizedSortedPValue.begin()).second,sequences,IMaps);
			BuildFrequencyMatrix(freqMatrixB, finderB, (*sequences.generalizedSortedPValue.begin()).second,background,IMaps);
			std::cout<<(*sequences.generalizedSortedPValue.begin()).second<<std::endl;
			bool foreground = true;
			PrintMap(freqMatrix,foreground);

			foreground=false;
			PrintMap(freqMatrixB,foreground);

			
			sequences.pValue=sequences.generalizedSortedPValue.begin()->first;
		}
		else{

			output.open("output.fasta",std::ios::out|std::ios::app);
			write(output,(*sequences.SortedPValue.begin()).second);
			write(output,"   ");
		
			output<<(*sequences.SortedPValue.begin()).first;
			write(output,"   ");
			output.close();
			BuildFrequencyMatrix(freqMatrix, finder, (*sequences.SortedPValue.begin()).second,sequences,IMaps);
			BuildFrequencyMatrix(freqMatrixB, finderB, (*sequences.SortedPValue.begin()).second,background,IMaps);
			std::cout<<(*sequences.SortedPValue.begin()).second<<std::endl;
			bool foreground = true;
			PrintMap(freqMatrix,foreground);

			foreground=false;
			PrintMap(freqMatrixB,foreground);

			
			sequences.pValue=sequences.SortedPValue.begin()->first;

		}

		/***
		PSSM
		***/
		BuildWeightMatrix(sequences,freqMatrix);
	
		
		freqMatrix.clear();
		freqMatrixB.clear();
		//std::cout<<"first "<<(*sequences.generalizedSortedPValue.begin()).first<<std::endl;
		/*****
			- Computes the probability of each nucleotide to appear at each position (from the top motif)
			- first output
		*****/
		//String<unsigned int> replaceString;
		//String<unsigned int> replaceStringB;
		
		

	


		

		sequences.generalizedKmer.clear();
		sequences.seqCounter.clear();
		sequences.SortedPValueReversed.clear();
		background.generalizedKmer.clear();
		background.seqCounter.clear();
		background.SortedPValue.clear();
		background.SortedPValueReversed.clear();
		
		//std::cout<<replaceString[0]<<" "<<replaceString[1]<<" "<<replaceString[2]<<" "<<leftBoundary(sequences.intervals[replaceString[0]][0])<<" "<<rightBoundary(sequences.intervals[replaceString[0]][0]);
		//replaceKmer(sequences,replaceString);
		//replaceKmer(background,replaceStringB);
		//clear(replaceString);
		//clear(replaceStringB);

		//PrintFastA(sequences);
		//clear(sequences.SArray);
		//clear(background.SArray);
		//std::cout<<sequences.generalizedSortedPValue.begin()->first<<std::endl;
		std::cout<<"vor if"<<std::endl;
		if(sequences.c==1){//nur im ersten Schritt creatIntervalTree danach addInterval

			for(unsigned s=0;s<sequences.SeqsNumber;++s){
				createIntervalTree(sequences.intervalTrees[s], sequences.intervals[s]);
			}
			for(unsigned b=0;b<background.SeqsNumber;++b){
				createIntervalTree(background.intervalTrees[b], background.intervals[b]);
			}
		}
		std::cout<<"nach if"<<std::endl;
		clear(sequences.intervals);
		clear(background.intervals);

		++sequences.c;
		++background.c;
	}
	while(sequences.pValue<0.05 && sequences.c<4);


	return 0;

}

