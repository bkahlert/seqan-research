
#include "dreme.h"

using namespace seqan;



int main(int argc, char const ** argv){
	
	Seq sequences;
    Seq background;
	IupacMaps IMaps;
	MapIupac(IMaps);//IupacMap for generalization
	char setting;


	if(argc <2 || argc>3){
		std::cerr<<"ERROR: Invalid argument count."<<std:: endl
				 <<"Usage:" <<argv[0]<<"File"<<std::endl;
		return 1;
	}


	
	readFastA(sequences,argv[1]);
	sequences.SeqsNumber=length(sequences.seqs);//number of sequences
	if(argc == 3)	readFastA(background,argv[2]);
	else{
		background.seqs=sequences.seqs;
		const int SEED = 7;
		Rng<MersenneTwister> rng(SEED);
		for(unsigned int i=0;i<sequences.SeqsNumber-1;++i){
	
			shuffle(background.seqs[i],rng);
			shuffle(background.seqs[i],rng);
			shuffle(background.seqs[i],rng);
			shuffle(background.seqs[i],rng);
		}
	
	}
	
	background.SeqsNumber=length(background.seqs);//number of sequences


	std::cout<<"-------------DREME-Settings-------------"<<std::endl;
	std::cout<<"----------------------------------------"<<std::endl;
	std::cout<<"Standard-Settings are: mink(6), maxk(6), Seed(100), pValueTreshold(0.05), NumberOfMotifs(5)"<<std::endl;
	while(true){
		std::cout<<"Alter the Settings? (Y) (N)"<<std::endl;
		std::cin>>setting;
		if(setting =='Y' || setting=='N' ||setting =='y' ||setting =='n') break;
		std::cout<<"**** ERROR: MUST BE Y OR N ****"<<std::endl;
	}
	unsigned kmer_len;
	unsigned kmer_len_end;
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
	
	while(true){
		std::cout<<"Save the data (PWMs) to a file? (Y) (N)"<<std::endl;
		std::cin>>setting;
		if(setting =='Y' || setting=='N' ||setting =='y' ||setting =='n') break;
		std::cout<<"**** ERROR: MUST BE Y OR N ****"<<std::endl;
	}
	bool save=false;
	if(setting=='Y' || setting =='y'){
		save=true;
		String<char> pwm;
		std::cout<<"Enter in the path, where the data should be saved:"<<std::endl;
		std::cin>>pwm;
		std::cout<<pwm;
	}
	

	/*std::ofstream outfile;
	outfile.open("test2.fasta");
	write2(outfile,sequences.ids,background.seqs,Fasta());

	outfile.close();*/
	
	std::ofstream output;
	std::ofstream PWM;
	
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
	priorFreq(background);

	
	
	/*****
			Start the Algorithm
	*****/
	do{
		std::cout<<"New Round "<<std::endl;
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
		InitGeneralization(IMaps,sequences,background);
		std::cout<<"Generalize done"<<std::endl;
		sequences.SortedPValueReversed.clear();

		
		
		
		Finder<Index<StringSet<String<Dna5> > > > finder(sequences.SArray);
		Finder<Index<StringSet<String<Dna5> > > > finderB(background.SArray);
		sequences.pValue=0;
		/*****
			- if there is not a single pValue<treshold (in the generalized motif list) compute the PWM of the motif with the smallest pValue
			--> else
		*****/
		if(sequences.generalizedSortedPValue.size()>0){
			
			std::map<String<Iupac>,unsigned int > seqCounter;
			std::map<String<Iupac>,unsigned int > backCounter;
			
	
			/*****
				- gets the top 100(generalized) motifs and computes the exact counter and pValue
			*****/
			exactGeneralizeCount(seqCounter,backCounter, finder, finderB,sequences,background, IMaps);
			std::cout<<"exactGeneralize done"<<std::endl;

			/****
					writes the motif with p-value into "output.fasta"
			****/
			output.open("output.fasta",std::ios::out|std::ios::app);
			write(output,(*sequences.generalizedSortedPValue.begin()).second);
			write(output,"   ");
		
			output<<(*sequences.generalizedSortedPValue.begin()).first;
			write(output,"   ");
			output.close();
			/****
				computes the PWM
			****/
			BuildFrequencyMatrix(finder, (*sequences.generalizedSortedPValue.begin()).second,sequences,IMaps);
			BuildFrequencyMatrix(finderB, (*sequences.generalizedSortedPValue.begin()).second,background,IMaps);
			std::cout<<(*sequences.generalizedSortedPValue.begin()).second<<std::endl;
			
			/****
					saves the PWM to a file
			****/
			saveData(sequences,PWM,sequences.c);
			
			sequences.pValue=sequences.generalizedSortedPValue.begin()->first;
		}
		else{
			/****
					writes the motif with p-value into "output.fasta"
			****/
			output.open("output.fasta",std::ios::out|std::ios::app);
			write(output,(*sequences.SortedPValue.begin()).second);
			write(output,"   ");
		
			output<<(*sequences.SortedPValue.begin()).first;
			write(output,"   ");
			output.close();
			/****
				computes the PWM
			****/
			BuildFrequencyMatrix( finder, (*sequences.SortedPValue.begin()).second,sequences,IMaps);
			BuildFrequencyMatrix( finderB, (*sequences.SortedPValue.begin()).second,background,IMaps);
			std::cout<<(*sequences.SortedPValue.begin()).second<<std::endl;
			/****
					saves the PWM to a file
			****/
			saveData(sequences,PWM,sequences.c);
			
			sequences.pValue=sequences.SortedPValue.begin()->first;

		}

		
		//BuildWeightMatrix(sequences);

		/***
		Information Content
		***/
		//BuildInformationContentMatrix(sequences);
	

		appendValue(sequences.allPWMs,sequences.freqMatrix);

		sequences.freqMatrix.clear();
		sequences.weightMatrix.clear();
		background.weightMatrix.clear();
		background.freqMatrix.clear();
		sequences.generalizedKmer.clear();
		sequences.seqCounter.clear();
		sequences.SortedPValueReversed.clear();
		background.generalizedKmer.clear();
		background.seqCounter.clear();
		background.SortedPValue.clear();
		background.SortedPValueReversed.clear();
		
	
		if(sequences.c==1){//if c==1 creatIntervalTree, otherwise addInterval

			for(unsigned s=0;s<sequences.SeqsNumber;++s){
				createIntervalTree(sequences.intervalTrees[s], sequences.intervals[s]);
			}
			for(unsigned b=0;b<background.SeqsNumber;++b){
				createIntervalTree(background.intervalTrees[b], background.intervals[b]);
			}
		}
		
		seqan::clear(sequences.intervals);
		seqan::clear(background.intervals);
		
		std::cout<<std::endl;
		++sequences.c;
		++background.c;
	}
	while(sequences.pValue<sequences.treshold && sequences.c<sequences.NumberOfMotifs);
	/******
			Compute the Clustering with GapStat
	******/
	String<Cluster> cluster;
	String<int> traceback;
	PWMClustering(sequences.allPWMs,cluster,traceback);
	std::cout<<"ComputeGapStat"<<std::endl<<std::endl;
	clusterTraceback(sequences.allMotifs,traceback,computeGapStat(sequences.allPWMs,cluster));
	
	seqan::clear(cluster);
	seqan::clear(traceback);
	seqan::clear(sequences.allMotifs);
	seqan::clear(sequences.allPWMs);
	
	return 0;

}

